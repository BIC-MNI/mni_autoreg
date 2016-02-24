/* ----------------------------- MNI Header -----------------------------------
@NAME       : init_params.c
@DESCRIPTION: collection of routines that will calculate the parameters necessary
              from an input transformation matrix for optimization when 
              mapping world coordinates  of volume 1 into world coords in volume 2.
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

@CREATED    : Thu May 27 16:50:50 EST 1993
                  
@MODIFIED   :  $Log: init_params.c,v $
@MODIFIED   :  Revision 96.12  2006-11-30 09:07:32  rotor
@MODIFIED   :   * many more changes for clean minc 2.0 build
@MODIFIED   :
@MODIFIED   :  Revision 96.11  2006/11/29 09:09:33  rotor
@MODIFIED   :   * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   :  Revision 96.10  2005/07/20 20:45:49  rotor
@MODIFIED   :      * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :      * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :      * Removed old VOLUME_IO cruft #defines
@MODIFIED   :      * Fixed up all Makefile.am's in subdirs
@MODIFIED   :      * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :      * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   :  Revision 96.9  2004/03/12 05:16:09  rotor
@MODIFIED   :   * fixed COV bug in init_params (lenezet)
@MODIFIED   :   * added epm-header.in for binary packaging
@MODIFIED   :   * Fixed -dircos problem in test cases
@MODIFIED   :
@MODIFIED   :  Revision 96.8  2004/02/12 05:54:27  rotor
@MODIFIED   :   * removed /static defs
@MODIFIED   :
@MODIFIED   :  Revision 96.7  2004/01/27 00:28:03  lenezet
@MODIFIED   :  change init_params to correct the COG bug when there is not input transform.
@MODIFIED   :  add the cosines director to the resampled field
@MODIFIED   :
@MODIFIED   :  Revision 96.6  2003/02/05 21:27:26  lenezet
@MODIFIED   :  array size correction.
@MODIFIED   :
@MODIFIED   :  Revision 96.5  2003/02/04 06:08:45  stever
@MODIFIED   :  Add support for correlation coefficient and sum-of-squared difference.
@MODIFIED   :
@MODIFIED   :  Revision 96.4  2002/12/13 21:16:30  lenezet
@MODIFIED   :  nonlinear in 2D has changed. The option -2D-non-lin is no more necessary. The grid transform has been adapted to feet on the target volume whatever is size. The Optimization is done on the dimensions for which "count" is greater than 1.
@MODIFIED   :
@MODIFIED   :  Revision 96.3  2002/08/14 19:54:26  lenezet
@MODIFIED   :   quaternion option added for the rotation
@MODIFIED   :
@MODIFIED   :  Revision 96.2  2002/03/26 14:15:40  stever
@MODIFIED   :  Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   :  Revision 96.1  1997/11/03 19:59:49  louis
@MODIFIED   :  - now include volume_io/internal_volume_io.h instead of volume_io.h
for(=; <; ++)
@MODIFIED   :
 * Revision 96.0  1996/08/21  18:21:58  louis
 * Release of MNI_AutoReg version 0.96
 *
 * Revision 9.6  1996/08/21  18:21:53  louis
 * Pre-release
 *
 * Revision 9.5  1996/08/12  14:15:50  louis
 * Never released version 0.95
 *
 * Revision 1.15  1996/08/12  14:15:44  louis
 * Pre-release
 *
 * Revision 1.14  1996/02/22  09:47:08  collins
 * Fixed loop limits for i and j when norming the principal axes (line 613)
 *
 * Revision 1.13  1995/09/11  12:37:16  collins
 * changes to init_transformation(): I not ensure that the principal
 * axes returned from cov_to_praxes form a righ-handed coordinate sytem.
 *
 * All refs to numerical recipes routines have been replaced.
 *
 * this is an updated working version - corresponds to mni_reg-0.1g
 *
 * Revision 1.12  1995/03/17  10:59:35  collins
 * corrected a memory bug - I was freeing two centroid vectors that were
 * not supposed to be freed in init_transformation.  These variables are
 * freed in the calling procedure init_params.
 *
 * Revision 1.11  1995/02/22  08:56:06  collins
 * Montreal Neurological Institute version.
 * compiled and working on SGI.  this is before any changes for SPARC/
 * Solaris.
 *
 * Revision 1.10  94/06/02  20:16:01  louis
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
 * Revision 1.9  94/04/26  12:54:19  louis
 * updated with new versions of make_rots, extract2_parameters_from_matrix 
 * that include proper interpretation of skew.
 * 
 * Revision 1.8  94/04/06  11:48:37  louis
 * working linted version of linear + non-linear registration based on Lvv
 * operator working in 3D
 * 
 * Revision 1.7  94/02/21  16:35:36  louis
 * version before feb 22 changes
 * 
 * Revision 1.6  93/11/15  16:26:44  louis
 * working version, with new library, with RCS revision stuff,
 * before deformations included
 * 
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /static-cvsroot/registration/mni_autoreg/minctracc/Numerical/init_params.c,v 96.12 2006-11-30 09:07:32 rotor Exp $";
#endif


#include "config.h"
#include <float.h>
#include <volume_io.h>

#include "constants.h"
#include "minctracc_arg_data.h"

#include "matrix_basics.h"
#include "cov_to_praxes.h"
#include "make_rots.h"
#include "quaternion.h"

extern Arg_Data *main_args;

#include "local_macros.h"
#include <Proglib.h>

#define  RAD_TO_DEG   (180.0 / M_PI)

VIO_BOOL rotmat_to_ang(float **rot, float *ang);

int point_not_masked(VIO_Volume volume, 
                            VIO_Real wx, VIO_Real wy, VIO_Real wz);

void set_up_lattice(VIO_Volume data,       /* in: volume  */
                           double *user_step, /* in: user requested spacing for lattice */
                           double *start,     /* out:world starting position of lattice  in volume dircos coords*/
                           double *wstart,     /* out:world starting position of lattice */
                           int    *count,     /* out:number of steps in each direction */
                           double *step,      /* out:step size in each direction */
                           VectorR directions[]);/* out: vector directions for each index*/


/* ----------------------------- MNI Header -----------------------------------
@NAME       : vol_cog - get the center of gravity of a volume.
@INPUT      : d1: one volume of data (already in memory).
              m1: its corresponding mask volume
              step: an 3 element array of step sizes in x,ya nd z directions
                
@OUTPUT     : centroid - vector giving centroid of points. This vector
                         must be defined by the calling routine.
@RETURNS    : TRUE if ok, FALSE if error.
@DESCRIPTION: this routine does calculates the cog using 
              volumetric subsampling in world space.
              These world coordinates are mapped back into each 
              data volume to get the actual value at the given location.

@GLOBALS    : 
@CALLS      : 
@CREATED    : Wed Aug  2 12:05:39 MET DST 1995 LC
@MODIFIED   : 
              
---------------------------------------------------------------------------- */
VIO_BOOL vol_cog(VIO_Volume d1, VIO_Volume m1, float *centroid, double *step)
{


  VectorR
    vector_step;

  PointR 
    starting_position,
    slice,
    row,
    col,
    voxel;

  VIO_Real
    tx,ty,tz;

  int
    i,
    r,c,s;

  float
    sx,sy,sz,si; 

  VIO_Real
    true_value;

  int
    count[ VIO_MAX_DIMENSIONS];
  double 
    start[VIO_MAX_DIMENSIONS],
    wstart[ VIO_MAX_DIMENSIONS],
    local_step[ VIO_MAX_DIMENSIONS];
  VectorR
    scaled_directions[VIO_MAX_DIMENSIONS],
    directions[ VIO_MAX_DIMENSIONS];  

                                /* build default sampling lattice info
                                   on the data set (d1)               */
  set_up_lattice(d1, step,start, wstart, count, local_step, directions);

  
   for(i=0; i<3; i++) {
     Point_x(scaled_directions[i]) = Point_x(directions[i]) * local_step[i];
     Point_y(scaled_directions[i]) = Point_y(directions[i]) * local_step[i];
     Point_z(scaled_directions[i]) = Point_z(directions[i]) * local_step[i];
   }

  fill_Point( starting_position, wstart[0], wstart[1], wstart[2]);
  
                                /* calculate centroids */

  sx = 0.0;
  sy = 0.0;
  sz = 0.0;
  si = 0.0;

 

  for(s=0; s<count[SLICE_IND]; s++) {

    SCALE_VECTOR( vector_step, scaled_directions[SLICE_IND], s);
    ADD_POINT_VECTOR( slice, starting_position, vector_step );

    for(r=0; r<count[ROW_IND]; r++) {

      SCALE_VECTOR( vector_step, scaled_directions[ROW_IND], r);
      ADD_POINT_VECTOR( row, slice, vector_step );

      SCALE_POINT( col, row, 1.0); /* init first col position */
      for(c=0; c<count[COL_IND]; c++) {

        convert_3D_world_to_voxel(d1, Point_x(col), Point_y(col), Point_z(col), &tx, &ty, &tz);


        fill_Point( voxel, tx, ty, tz ); /* build the voxel POINT */
        

        if (point_not_masked(m1, Point_x(col), Point_y(col), Point_z(col))) {        
          
          if (INTERPOLATE_TRUE_VALUE( d1, &voxel, &true_value )) {
            

            sx +=  Point_x(col) * true_value;
            sy +=  Point_y(col) * true_value;
            sz +=  Point_z(col) * true_value;
            
            si += true_value;
          }
          /* else requested voxel is just outside volume., so ignore it */

        }
        
        ADD_POINT_VECTOR( col, col, scaled_directions[COL_IND] );
        
      }
    }
  }

  if (si!=0.0) {
    centroid[1] = sx/ si;
    centroid[2] = sy/ si;
    centroid[3] = sz/ si;
    
    return(TRUE);
    
  }
  else {
    return(FALSE);
  }
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : vol_cov - get the covariance of a data volume
@INPUT      : d1: one volume of data (already in memory).
              m1: its corresponding mask volume
              step: an 3 element array of step sizes in x,ya nd z directions
              centroid - vector giving centroid of points. This vector
                         must be defined by the calling routine.  
@OUTPUT     : covar    - covariance matrix (in zero offset form).
                         must be defined by the calling routine.
@RETURNS    : TRUE if ok, FALSE if error.
@DESCRIPTION: this routine does calculates the covariance using 
              volumetric subsampling in world space.
              These world coordinates are mapped back into each 
              data volume to get the actual value at the given location.
              I use the header info and kernel size to calculate the
              positions of the sub-samples in each vol.

@GLOBALS    : 
@CALLS      : 
@CREATED    : Wed Aug  2 12:05:39 MET DST 1995 LC
@MODIFIED   : 
              
---------------------------------------------------------------------------- */
VIO_BOOL vol_cov(VIO_Volume d1, VIO_Volume m1, float *centroid, float **covar, double *step)
{


  VectorR
    vector_step;

  PointR 
    starting_position,
    slice,
    row,
    col,
    voxel;

  VIO_Real
    tx,ty,tz;

  int
    i,r,c,s;

  float
    sxx,syy,szz,
    sxy,syz,sxz,
    si; 

  VIO_Real
    true_value;

  int
    count[VIO_MAX_DIMENSIONS];
  double 
    start[VIO_MAX_DIMENSIONS],
    wstart[VIO_MAX_DIMENSIONS],
    local_step[VIO_MAX_DIMENSIONS];
  VectorR
    scaled_directions[VIO_MAX_DIMENSIONS],
    directions[VIO_MAX_DIMENSIONS];  

                                /* build default sampling lattice info
                                   on the data set (d1)               */

  set_up_lattice(d1, step,
                 start, wstart, count, local_step, directions);

   for(i=0; i<3; i++) {
     Point_x(scaled_directions[i]) = Point_x(directions[i]) * local_step[i];
     Point_y(scaled_directions[i]) = Point_y(directions[i]) * local_step[i];
     Point_z(scaled_directions[i]) = Point_z(directions[i]) * local_step[i];
   }

  
  fill_Point( starting_position, wstart[0], wstart[1], wstart[2]);
  
  si = 0.0;
  sxx = syy = szz = 0.0;
  sxy = syz = sxz = 0.0;
    
                                /* now calculate variances and co-variances */

  for(s=0; s<count[SLICE_IND]; s++) {
    
    SCALE_VECTOR( vector_step, scaled_directions[SLICE_IND], s);
    ADD_POINT_VECTOR( slice, starting_position, vector_step );
    
    for(r=0; r<count[ROW_IND]; r++) {
      
      SCALE_VECTOR( vector_step, scaled_directions[ROW_IND], r);
      ADD_POINT_VECTOR( row, slice, vector_step );
      
      SCALE_POINT( col, row, 1.0); /* init first col position */
      for(c=0; c<count[COL_IND]; c++) {
        
        
        convert_3D_world_to_voxel(d1, Point_x(col), Point_y(col), Point_z(col), &tx, &ty, &tz);
        
        fill_Point( voxel, tx, ty, tz ); /* build the voxel POINT */
        
        if (point_not_masked(m1, Point_x(col), Point_y(col), Point_z(col))) {        
          
          if (INTERPOLATE_TRUE_VALUE( d1, &voxel, &true_value )) {
            
            sxx += (Point_x( col )-centroid[1]) * (Point_x( col )-centroid[1]) * true_value;
            syy += (Point_y( col )-centroid[2]) * (Point_y( col )-centroid[2]) * true_value;
            szz += (Point_z( col )-centroid[3]) * (Point_z( col )-centroid[3]) * true_value;
            sxy += (Point_x( col )-centroid[1]) * (Point_y( col )-centroid[2]) * true_value;
            syz += (Point_y( col )-centroid[2]) * (Point_z( col )-centroid[3]) * true_value;
            sxz += (Point_x( col )-centroid[1]) * (Point_z( col )-centroid[3]) * true_value;
            si += true_value;
          }
          /* else requested voxel is just outside volume., so ignore it */
          
        } 
        
        ADD_POINT_VECTOR( col, col, scaled_directions[COL_IND] );
        
      }
    }
  }

  if (si != 0.0) {
    covar[1][1] = sxx/si; covar[1][2] = sxy/si; covar[1][3] = sxz/si;
    covar[2][1] = sxy/si; covar[2][2] = syy/si; covar[2][3] = syz/si;
    covar[3][1] = sxz/si; covar[3][2] = syz/si; covar[3][3] = szz/si;
    
    return(TRUE);
    
  }
  else {
    return(FALSE);
  }
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : vol_to_cov - get covariance and cog of volume.
@INPUT      : d1: one volume of data (already in memory).
              m1: its corresponding mask volume
              step: an 3 element array of step sizes in x,ya nd z directions
                
@OUTPUT     : centroid - vector giving centroid of points. This vector
                         must be defined by the calling routine.
              covar    - covariance matrix (in zero offset form).
                         must be defined by the calling routine.
@RETURNS    : TRUE if ok, FALSE if error.
@DESCRIPTION: this routine does calculates the covariance using 
              volumetric subsampling in world space.
              These world coordinates are mapped back into each 
              data volume to get the actual value at the given location.
              I use the header info and kernel size to calculate the
              positions of the sub-samples in each vol.

@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb 5, 1992 lc
@MODIFIED   : Thu May 27 16:50:50 EST 1993 lc
                 rewrite for minc files and david's library
---------------------------------------------------------------------------- */
VIO_BOOL vol_to_cov(VIO_Volume d1, VIO_Volume m1, float *centroid, float **covar, double *step)
{
  int
    i,count[VIO_MAX_DIMENSIONS];
  double 
    start[VIO_MAX_DIMENSIONS],
    wstart[VIO_MAX_DIMENSIONS],
    local_step[VIO_MAX_DIMENSIONS];
  VectorR
    directions[VIO_MAX_DIMENSIONS];  

  if (main_args->flags.debug) {

    set_up_lattice(d1, step,
                   start, wstart, count, local_step, directions);

    print ("in vol to cov\n");
    print ("start = %8.2f %8.2f %8.2f \n",start[0],start[1],start[2]);
    print ("count = %8d %8d %8d \n",count[0],count[1],count[2]);
    print ("step  = %8.2f %8.2f %8.2f \n",local_step[0],local_step[1],local_step[2]);
    
    for(i=0; i<3; i++)
      print ("direct= %8.2f %8.2f %8.2f \n",
             Point_x(directions[i]),
             Point_y(directions[i]),
             Point_z(directions[i]));
  }


  if ( vol_cog(d1, m1, centroid, step) )
    
    return ( vol_cov( d1, m1, centroid, covar, step) );

  else
    
    return (FALSE);

}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : init_transformation - get trans parameters using principal axis
                 transformation technique.  the transformation points in 
                 d1 -> d2.
@INPUT      : d1,d2:
                two volumes of data (already in memory).
              m1,m2:
                two mask volumes for data (already in memory).
              step: an 3 element array of step sizes in x,y and z directions
              verbose: =0 for quiet, >1 otherwise
                
@OUTPUT     : c1,c2    - vector giving volume centroids. This vector
                         must be defined by the calling routine.

                         c1 will not be overwritten if it already contains valid data.

              trans    - translation matrix (in zero offset form).
                         must be defined by the calling routine.
              rots     - rotation matrix (in zero offset form).
                         must be defined by the calling routine.
              ang      - vector of rotation angles, in radians (in zero offset form).
                         must be defined by the calling routine.
              c1       - vector for centroid of d1 (in zero offset form).
                         must be defined by the calling routine.
              c2       - vector for centroid of d2 (in zero offset form).
                         must be defined by the calling routine.
@RETURNS    : TRUE if ok, FALSE if error.
@DESCRIPTION: this routine uses the principal axis method to determine the 
              initial guess for registration transformation
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb 5, 1992 lc
@MODIFIED   : Thu May 27 16:50:50 EST 1993 lc
                 rewrite for minc files and david's library
Wed Aug  2 12:05:39 MET DST 1995 LC
  c1 is not overwritten if it contains valid data.  ie I won't recalculate the
  COG of volume 1 if the values in C1 are not==-DBL_MAX.
---------------------------------------------------------------------------- */
static  VIO_BOOL init_transformation(
                                     VIO_Volume d1, /* data for volume1 */
                                     VIO_Volume d2, /* data for volume2 */
                                     VIO_Volume m1, /* mask for volume1 */
                                     VIO_Volume m2, /* mask for volume2 */
                                     double *step, /* in x,y,z order  */
                                     int    verbose,
                                     float **trans,     /* translation matrix to go from d1 to d2 */
                                     float **rots,      /* rotation matrix to go from d1 to d2    */
                                     float *ang,        /* rotation angles to go from d1 to d2    */
                                     float *c1,         /* centroid of masked d1 */
                                     float *c2,         /* centroid of masked d1 */
                                     float *scale,      /* scaling from d1 to d2 */
                                     int forced_center,
                                     Transform_Flags *flags) /* flags for estimation */
{
  float
    dir,
    tx,ty,tz,
    *angles,                        /* rotation angles - rx,ry and rz */
    **cov1,**cov2,                /* covariance matrix */
    **prin_axes1, **prin_axes2, /* principal axis */
    **R1,**R2,                        /* rotation matrix (normalized prin_axes) */
    **Rinv,**R;
  
  float
    norm;
  
  int
    stat,
    ndim,i,j;

  nr_identf(trans,1,4,1,4);        /* start with identity                       */
  nr_identf(rots, 1,4,1,4);
  
  VIO_ALLOC2D(cov1       ,4,4);
  VIO_ALLOC2D(cov2       ,4,4);
  VIO_ALLOC2D(prin_axes1 ,4,4);
  VIO_ALLOC2D(prin_axes2 ,4,4);
  VIO_ALLOC2D(R1         ,4,4);
  VIO_ALLOC2D(R2         ,4,4);
  VIO_ALLOC2D(R          ,4,4);
  VIO_ALLOC2D(Rinv       ,4,4);
  ALLOC(angles     ,4);
  

  stat = TRUE;

  /* =========  calculate COG and COV for volume 1   =======  */
                                /* if center already set, then don't recalculate */
  if ( !forced_center) {
    stat = vol_cog(d1, m1, c1, step);
    if (verbose>0 && stat) print ("COG of v1: %f %f %f\n",c1[1],c1[2],c1[3]);
  }
  else {
    if (verbose>0) print ("COG of v1 forced: %f %f %f\n",c1[1],c1[2],c1[3]);
  }
  if (!stat || !vol_cov(d1, m1, c1, cov1, step ) ) {
    print_error_and_line_num("%s", __FILE__, __LINE__,"Cannot calculate the COG or COV of volume 1.\n" );
    return(FALSE);
  }


  /* =========  calculate COG and COV for volume 2 only if needed:   =======  */

  if (flags->estimate_trans || flags->estimate_rots || flags->estimate_scale) {
    if (! vol_to_cov(d2, m2, c2, cov2, step ) ) {
      print_error_and_line_num("%s", __FILE__, __LINE__,"Cannot calculate the COG or COV of volume 2.\n" );
      return(FALSE);
    }
    if (verbose>0) print ("COG of v2: %f %f %f\n",c2[1],c2[2],c2[3]);
  }
  else {
    if (verbose>0) print ("Only center required, now returning from init_transformation\n");
  }

  if (flags->estimate_trans) {
    tx = c2[1] - c1[1];    /* translations to map vol1 into vol2                  */
    ty = c2[2] - c1[2];
    tz = c2[3] - c1[3];

    if (verbose>0) print ("   [trans] = %f %f %f\n",tx,ty,tz);

  }
  else {
    tx = ty = tz = 0.0;
  }


  trans[1][4] += tx;           /* set translations in translation matrix        */
  trans[2][4] += ty;
  trans[3][4] += tz;
  
  if (flags->estimate_rots || flags->estimate_scale) {

                                /* get the principal axes, returned in
                                   cols of prin_axes{1,2} */

    cov_to_praxes(3, cov1, prin_axes1);   
    cov_to_praxes(3, cov2, prin_axes2);

    if (verbose > 1) {
      print ("cov1:                              cov2:\n");
      for (i=1; i<=3; i++) {
        for (j=1; j<=3; j++)
          print ("%8.3f ", cov1[i][j]);
        print ("|");
        for (j=1; j<=3; j++)
          print ("%8.3f ", cov2[i][j]);
        print ("\n\n");
      }
    }
    
                                /* make sure that both sets of
                                   principal axes represent
                                   right-handed coordinate system:
                                   [(p1 x p1).p3]>0, where 'x' is
                                   cross product and '.' is dot
                                   product*/

    dir = prin_axes1[1][3] * (prin_axes1[2][1] * prin_axes1[3][2] - 
                              prin_axes1[3][1] * prin_axes1[2][2])  +
          prin_axes1[2][3] * (prin_axes1[1][1] * prin_axes1[3][2] - 
                              prin_axes1[3][1] * prin_axes1[1][2])  +
          prin_axes1[3][3] * (prin_axes1[1][1] * prin_axes1[2][2] - 
                              prin_axes1[2][1] * prin_axes1[1][2]);
    if (dir < 0) {                /* if lefthanded, change dir of 3rd vector */
      prin_axes1[2][3] *= -1.0;
      prin_axes1[3][3] *= -1.0;
    }

    dir = prin_axes2[1][3] * (prin_axes2[2][1] * prin_axes2[3][2] - 
                              prin_axes2[3][1] * prin_axes2[2][2])  +
          prin_axes2[2][3] * (prin_axes2[1][1] * prin_axes2[3][2] - 
                              prin_axes2[3][1] * prin_axes2[1][2])  +
          prin_axes2[3][3] * (prin_axes2[1][1] * prin_axes2[2][2] - 
                              prin_axes2[2][1] * prin_axes2[1][2]);
    if (dir < 0) {                /* if lefthanded, change dir of 3rd vector */
      prin_axes2[2][3] *= -1.0;
      prin_axes2[3][3] *= -1.0;
    }

                                /* print out the prin axes */
    if (verbose > 1) {
      print ("prin_axes1:                      princ_axes2:\n");
      for (i=1; i<=3; i++) {
        for (j=1; j<=3; j++)
          print ("%8.3f ", prin_axes1[i][j]);
        print ("|");
        for (j=1; j<=3; j++)
          print ("%8.3f ", prin_axes2[i][j]);
        print ("\n\n");
      }
    }

                                /* build rotation matrixes from principal axes */

    nr_copyf(prin_axes1,1,3,1,3,R1);
    nr_copyf(prin_axes2,1,3,1,3,R2);
    
    
                                /* normalize the principal axes lengths */
    for (j=1; j<=3; ++j) {        
      norm = sqrt( prin_axes1[1][j]*prin_axes1[1][j] + 
                   prin_axes1[2][j]*prin_axes1[2][j] + 
                   prin_axes1[3][j]*prin_axes1[3][j]);
      for (i=1; i<=3; ++i)
        R1[i][j] /= norm;
      
      norm = sqrt( prin_axes2[1][j]*prin_axes2[1][j] + 
                   prin_axes2[2][j]*prin_axes2[2][j] + 
                   prin_axes2[3][j]*prin_axes2[3][j]);
      for (i=1; i<=3; ++i)
        R2[i][j] /= norm;
    }
    
    invertmatrix(3,R1,Rinv);
    
    nr_multf(Rinv,1,3,1,3, R2,1,3,1,3, R);
    
    
    if (verbose > 1) {
      print ("r1:                         r2:                        r:\n");
      for (i=1; i<=3; i++) {
        for (j=1; j<=3; j++)
          print ("%8.3f ", R1[i][j]);
        print ("|");
        for (j=1; j<=3; j++)
          print ("%8.3f ", R2[i][j]);
        print ("|");
        for (j=1; j<=3; j++)
          print ("%8.3f ", R[i][j]);
        print ("\n");
      }
    }
    
    
    transpose(3,3,R,R);                /* all of the princ axes stuff uses vec*mat */
    
    if (!rotmat_to_ang(R, angles)) {
      (void)fprintf(stderr,"Could not extract angles from rotation matrix:\n");
      printmatrix(3,3,R);

      (void)fprintf(stderr,"You can try rerunning minctracc with '-est_translation' \n");
      (void)fprintf(stderr,"The program will continue with angles = 0.0, 0.0, 0.0   \n");
      angles[1] =angles[2] = angles[3] = 0.0;

    }

 

    scale[1] = 1.0;
    scale[2] = 1.0;
    scale[3] = 1.0;
    
    ang[1] = angles[1];                /* rotation about X axis                   */
    ang[2] = angles[2];                /* rotation about Y axis                   */
    ang[3] = angles[3];                /* rotation about Z axis                   */
    
    for (i=1; i<=3; ++i)        /* set rotations in matrix                 */
      for (j=1; j<=3; ++j) {
        rots[i][j] = R[i][j];
      }
   

    ndim = 3;
    
    if (verbose > 1) {
      (void) print("\nFor volume 1 :");
      (void) print("\nCentroid :");
      for (i=1; i<=ndim; i++) (void) print("  %f",c1[i]);
      (void) print("\n\n");
      (void) print("Principal axes\n");
      for (i=1; i<=ndim; i++) {
        (void) print("VIO_Vector %d :",i);
        for (j=1; j<=ndim; j++) {
          (void) print("  %f",prin_axes1[j][i]);
        }
        (void) print("\n");
      }
      
      (void) print("\n");
      
      (void) print("\nFor volume 2 :");
      (void) print("\nCentroid :");
      for (i=1; i<=ndim; i++) (void) print("  %f",c2[i]);
      (void) print("\n\n");
      (void) print("Principal axes\n");
      for (i=1; i<=ndim; i++) {
        (void) print("VIO_Vector %d :",i);
        for (j=1; j<=ndim; j++) {
          (void) print("  %f",prin_axes2[j][i]);
        }
        (void) print("\n");
      }
      
      (void) print ("rotation angles are: %f %f %f\n",
                    angles[1]*RAD_TO_DEG,
                    angles[2]*RAD_TO_DEG,
                    angles[3]*RAD_TO_DEG);
      (void) print ("translation mm     : %f %f %f\n",tx,ty,tz);
    }
    
  }
  else {

    if (verbose>0) print ("Only center & trans required, now returning from init_transformation\n");

    scale[1] = 1.0;
    scale[2] = 1.0;
    scale[3] = 1.0;
    
    ang[1] = 0.0;
    ang[2] = 0.0;
    ang[3] = 0.0;
  }

  VIO_FREE2D(cov1);
  VIO_FREE2D(cov2);
  VIO_FREE2D(prin_axes1);
  VIO_FREE2D(prin_axes2);
  VIO_FREE2D(R1);
  VIO_FREE2D(R2);
  VIO_FREE2D(R);
  VIO_FREE2D(Rinv);
  FREE(angles);
  
  return(TRUE);
}




/*********************************************************************************************************



                                   init_transformation_quater()

same as init_transformation but with quaternions

@OUTPUT     : instead of rots it is quats
              quats     - rotation matrix with quaternion
                         must be defined by the calling routine.

CREATED@ 1 May 2002

*********************************************************************************************************/

static  VIO_BOOL init_transformation_quater(
                                            VIO_Volume d1, /* data for volume1 */
                                            VIO_Volume d2, /* data for volume2 */
                                            VIO_Volume m1, /* mask for volume1 */
                                            VIO_Volume m2, /* mask for volume2 */
                                            double *step, /* in x,y,z order  */
                                            int    verbose,
                                            float **trans,     /* translation matrix to go from d1 to d2 */
                                            float *ang,        /* rotation angles to go from d1 to d2    */
                                            float *quats,      /* quaternions to go from d1 to d2 */
                                            float *c1,         /* centroid of masked d1 */
                                            float *c2,         /* centroid of masked d1 */
                                            float *scale,      /* scaling from d1 to d2 */
                                            int forced_center,
                                            Transform_Flags *flags) /* flags for estimation */
{
  float
    dir,
    tx,ty,tz,
    **cov1,**cov2,                /* covariance matrix */
    **prin_axes1, **prin_axes2, /* principal axis */
    **R1,**R2,                        /* rotation matrix (normalized prin_axes) */
    **Rinv,**R;
  
  double *qt,*vec,phi;
  
  float
    norm;
  
  int
    stat,
    ndim,i,j;

  
  
  VIO_ALLOC2D(cov1       ,4,4);
  VIO_ALLOC2D(cov2       ,4,4);
  VIO_ALLOC2D(prin_axes1 ,4,4);
  VIO_ALLOC2D(prin_axes2 ,4,4);
  VIO_ALLOC2D(R1         ,4,4);
  VIO_ALLOC2D(R2         ,4,4);
  VIO_ALLOC2D(R          ,4,4);
  VIO_ALLOC2D(Rinv       ,4,4);
  ALLOC(qt,4);
  
  nr_identf(trans,1,4,1,4);        /* start with identity                       */
  for(i=0; i<3; i++) qt[i]=0.0;
  qt[3]=1.0;

  stat = TRUE;

  /* =========  calculate COG and COV for volume 1   =======  */

                                /* if center already set, then don't recalculate */
  if ( !forced_center) {
    stat = vol_cog(d1, m1, c1, step);
    if (verbose>0 && stat) print ("COG of v1: %f %f %f\n",c1[1],c1[2],c1[3]);
  }
  else {
    if (verbose>0) print ("COG of v1 forced: %f %f %f\n",c1[1],c1[2],c1[3]);
  }

  if (!stat || !vol_cov(d1, m1, c1, cov1, step ) ) {
    print_error_and_line_num("%s", __FILE__, __LINE__,"Cannot calculate the COG or COV of volume 1\n." );
    return(FALSE);
  }

  /* =========  calculate COG and COV for volume 2 only if needed:   =======  */

  if (flags->estimate_trans || flags->estimate_quats || flags->estimate_scale) {
    if (! vol_to_cov(d2, m2, c2, cov2, step ) ) {
      print_error_and_line_num("%s", __FILE__, __LINE__,"Cannot calculate the COG or COV of volume 2\n." );
      return(FALSE);
    }
    if (verbose>0) print ("COG of v2: %f %f %f\n",c2[1],c2[2],c2[3]);
  }
  else {
    if (verbose>0) print ("Only center required, now returning from init_transformation\n");
  }

  if (flags->estimate_trans) {
    tx = c2[1] - c1[1];    /* translations to map vol1 into vol2                  */
    ty = c2[2] - c1[2];
    tz = c2[3] - c1[3];

    if (verbose>0) print ("   [trans] = %f %f %f\n",tx,ty,tz);

  }
  else {
    tx = ty = tz = 0.0;
  }


  trans[1][4] += tx;           /* set translations in translation matrix        */
  trans[2][4] += ty;
  trans[3][4] += tz;

  if (flags->estimate_quats || flags->estimate_scale) {

                                /* get the principal axes, returned in
                                   cols of prin_axes{1,2} */

    cov_to_praxes(3, cov1, prin_axes1);   
    cov_to_praxes(3, cov2, prin_axes2);

    if (verbose > 1) {
      print ("cov1:                              cov2:\n");
      for (i=1; i<=3; i++) {
        for (j=1; j<=3; j++)
          print ("%8.3f ", cov1[i][j]);
        print ("|");
        for (j=1; j<=3; j++)
          print ("%8.3f ", cov2[i][j]);
        print ("\n\n");
      }
    }
    
                                /* make sure that both sets of
                                   principal axes represent
                                   right-handed coordinate system:
                                   [(p1 x p1).p3]>0, where 'x' is
                                   cross product and '.' is dot
                                   product*/

    dir = prin_axes1[1][3] * (prin_axes1[2][1] * prin_axes1[3][2] - 
                              prin_axes1[3][1] * prin_axes1[2][2])  +
          prin_axes1[2][3] * (prin_axes1[1][1] * prin_axes1[3][2] - 
                              prin_axes1[3][1] * prin_axes1[1][2])  +
          prin_axes1[3][3] * (prin_axes1[1][1] * prin_axes1[2][2] - 
                              prin_axes1[2][1] * prin_axes1[1][2]);
    if (dir < 0) {                /* if lefthanded, change dir of 3rd vector */
      prin_axes1[2][3] *= -1.0;
      prin_axes1[3][3] *= -1.0;
    }

    dir = prin_axes2[1][3] * (prin_axes2[2][1] * prin_axes2[3][2] - 
                              prin_axes2[3][1] * prin_axes2[2][2])  +
          prin_axes2[2][3] * (prin_axes2[1][1] * prin_axes2[3][2] - 
                              prin_axes2[3][1] * prin_axes2[1][2])  +
          prin_axes2[3][3] * (prin_axes2[1][1] * prin_axes2[2][2] - 
                              prin_axes2[2][1] * prin_axes2[1][2]);
    if (dir < 0) {                /* if lefthanded, change dir of 3rd vector */
      prin_axes2[2][3] *= -1.0;
      prin_axes2[3][3] *= -1.0;
    }

                                /* print out the prin axes */
    if (verbose > 1) {
      print ("prin_axes1:                      princ_axes2:\n");
      for (i=1; i<=3; i++) {
        for (j=1; j<=3; j++)
          print ("%8.3f ", prin_axes1[i][j]);
        print ("|");
        for (j=1; j<=3; j++)
          print ("%8.3f ", prin_axes2[i][j]);
        print ("\n\n");
      }
    }

                                /* build rotation matrixes from principal axes */

    nr_copyf(prin_axes1,1,3,1,3,R1);
    nr_copyf(prin_axes2,1,3,1,3,R2);
    
    
                                /* normalize the principal axes lengths */
    for (j=1; j<=3; ++j) {        
      norm = sqrt( prin_axes1[1][j]*prin_axes1[1][j] + 
                   prin_axes1[2][j]*prin_axes1[2][j] + 
                   prin_axes1[3][j]*prin_axes1[3][j]);
      for (i=1; i<=3; ++i)
        R1[i][j] /= norm;
      
      norm = sqrt( prin_axes2[1][j]*prin_axes2[1][j] + 
                   prin_axes2[2][j]*prin_axes2[2][j] + 
                   prin_axes2[3][j]*prin_axes2[3][j]);
      for (i=1; i<=3; ++i)
        R2[i][j] /= norm;
    }
    
    invertmatrix(3,R1,Rinv);
    
    nr_multf(Rinv,1,3,1,3, R2,1,3,1,3, R);
    
    
    if (verbose > 1) {

      print ("r1:                         r2:                        r:\n");
      for (i=1; i<=3; i++) {
        for (j=1; j<=3; j++)
          print ("%8.3f ", R1[i][j]);
        print ("|");
        for (j=1; j<=3; j++)
          print ("%8.3f ", R2[i][j]);
        print ("|");
        for (j=1; j<=3; j++)
          print ("%8.3f ", R[i][j]);
        print ("\n");
      }
    }
    
    
    transpose(3,3,R,R);                /* all of the princ axes stuff uses vec*mat */

    /* extract quaternion from the matrix that we have find */
    extract_quaternions(R, qt);


    scale[1] = 1.0;
    scale[2] = 1.0;
    scale[3] = 1.0;
    
    ndim = 3;

    if (verbose > 1) {
      (void) print("\nFor volume 1 :");
      (void) print("\nCentroid :");
      for (i=1; i<=ndim; i++) (void) print("  %f",c1[i]);
      (void) print("\n\n");
      (void) print("Principal axes\n");
      for (i=1; i<=ndim; i++) {
        (void) print("VIO_Vector %d :",i);
        for (j=1; j<=ndim; j++) {
          (void) print("  %f",prin_axes1[j][i]);
        }
        (void) print("\n");
      }
      
      (void) print("\n");
      
      (void) print("\nFor volume 2 :");
      (void) print("\nCentroid :");
      for (i=1; i<=ndim; i++) (void) print("  %f",c2[i]);
      (void) print("\n\n");
      (void) print("Principal axes\n");
      for (i=1; i<=ndim; i++) {
        (void) print("VIO_Vector %d :",i);
        for (j=1; j<=ndim; j++) {
          (void) print("  %f",prin_axes2[j][i]);
        }
        (void) print("\n");
      }
      ALLOC(vec,3);
      quat_to_axis(vec,&phi,qt);
      
      (void) print ("\n\nrotation vector : %f %f %f the rotation angle is %f deg\n",
                    vec[1],
                    vec[2],
                    vec[3], phi);
      FREE(vec);
      (void) print ("translation mm     : %f %f %f\n",tx,ty,tz);
    }
    
  }
  else {

    if (verbose>0) print ("Only center & trans required, now returning from init_transformation\n");

    scale[1] = 1.0;
    scale[2] = 1.0;
    scale[3] = 1.0;
    
    ang[1] = 0.0;
    ang[2] = 0.0;
    ang[3] = 0.0;
  }
  
  for(i=0; i<3; i++) quats[i]=qt[i];

  VIO_FREE2D(cov1);
  VIO_FREE2D(cov2);
  VIO_FREE2D(prin_axes1);
  VIO_FREE2D(prin_axes2);
  VIO_FREE2D(R1);
  VIO_FREE2D(R2);
  VIO_FREE2D(R);
  VIO_FREE2D(Rinv);
  FREE(qt);
  
  return(TRUE);
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : init_params
                get the parameters necessary to map volume 1 to volume 2
@INPUT      : d1,d2:
                two volumes of data (already in memory).
              m1,m2:
                two mask volumes for data (already in memory).
              globals:
                a global data structure containing info from the command line.
                
@OUTPUT     : 
@RETURNS    : TRUE if ok, FALSE if error.
@DESCRIPTION: 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Thu May 27 16:50:50 EST 1993 lc
                 written for minc files and david's library
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_BOOL init_params(VIO_Volume d1,
                 VIO_Volume d2,
                 VIO_Volume m1,
                 VIO_Volume m2, 
                 Arg_Data *globals)
{  
  float 
    **trans,                        /* principal componant variables */
    **rots,  
    *quats,
    **cov1,
    *c1,
    *c2,
    *ang,
    *sc,
    quats3;
    
  int
    center_forced,i;
    
  VIO_Transform
    *lt;

  /* if there is no input transformation information specified on the
     command line, then the principal axis transformation (PAT) is
     used to initialize the parameters that will be used in the
     optimization process,

     otherwise,

     the command line transformation info must be interpreted to
     initialize the parameters that will be optimized.
  */


                                /* set a flag to see if the center was
                                   set on the command line */


  center_forced = (globals->trans_info.center[0] != -DBL_MAX || 
                   globals->trans_info.center[1] != -DBL_MAX || 
                   globals->trans_info.center[2] != -DBL_MAX);
  
  if (globals->flags.debug && center_forced) {
    print ("Center of rot/scale forced: %f %f %f\n",
           globals->trans_info.center[0],
           globals->trans_info.center[1],
           globals->trans_info.center[2]);
  }

                                /* if no transformation specified on the command line,
                                   of if PAT transformation selected,
                                     do principal axes transformation to init trans. */

  if (globals->trans_info.use_default == TRUE || 
      globals->trans_info.transform_type==TRANS_PAT) {

    if (globals->flags.debug) print ("  using PAT to get initial parameters:\n");
                                /* if 
                                      none of the -est_* flags are set on the command line,
                                      or if PAT transformation requested, 
                                   then 
                                       assume all -est_* to be true so that they will
                                       estimated in init_transformation and stored, 
                                   otherwise 
                                       leave the -est_* as set on the command line,
                                       and store only the ones requested. */

      if ((globals->trans_info.transform_type==TRANS_PAT) || 
          !(strlen(globals->filenames.measure_file)!=0 ||
            globals->trans_flags.estimate_center ||
            globals->trans_flags.estimate_scale  ||
            globals->trans_flags.estimate_rots   ||
            globals->trans_flags.estimate_trans)) {
        globals->trans_flags.estimate_center = TRUE;
        globals->trans_flags.estimate_scale = TRUE; 
        if(globals->trans_info.rotation_type == TRANS_ROT)
          {
            globals->trans_flags.estimate_rots = TRUE;
            globals->trans_flags.estimate_quats = FALSE;
          }
        else
          {
            globals->trans_flags.estimate_rots = FALSE;
            globals->trans_flags.estimate_quats = TRUE;
          }
        globals->trans_flags.estimate_trans = TRUE;
      }
      
      if (globals->flags.debug) 
        {
          print ("  will try to get:");
          if (globals->trans_flags.estimate_center) print (" [center]");
          if (globals->trans_flags.estimate_scale)  print (" [scale]");
          if (globals->trans_flags.estimate_rots)   print (" [rots]");
          if (globals->trans_flags.estimate_quats)   print (" [quats]");
          if (globals->trans_flags.estimate_trans)  print (" [trans]");
          print ("\n");
        }
      
      /* -est_* flags are now set, continue with the PAT */
      
      VIO_ALLOC2D(trans,5,5);
      VIO_ALLOC2D(rots,5,5);
      ALLOC(ang ,4);
      ALLOC(quats ,4);
      ALLOC(c1 ,4);
      ALLOC(c2 ,4);
      ALLOC(sc ,4);
    
                                /* set c1 to the value possibly forced on the command line */
      for(i=0; i<3; i++) c1[i+1] = globals->trans_info.center[i];
      if(globals->trans_info.rotation_type == TRANS_ROT)
        {
          if (!init_transformation(d1,d2,m1,m2, globals->step, globals->flags.verbose,
                                   trans,rots,ang,c1,c2,sc, center_forced,
                                   &(globals->trans_flags)))
            return(FALSE);
        }
      else
        {
          if (!init_transformation_quater(d1,d2,m1,m2, globals->step, globals->flags.verbose,
                                          trans,ang,quats,c1,c2,sc, center_forced,
                                          &(globals->trans_flags)))
            return(FALSE);
        }
      
      for(i=0; i<3; i++) {
        if (globals->trans_flags.estimate_rots)
          globals->trans_info.rotations[i]    = ang[i+1];
        else
          globals->trans_info.rotations[i]    = 0.0;
        
        if (globals->trans_flags.estimate_trans)
          globals->trans_info.translations[i] = trans[i+1][4];
        else
          globals->trans_info.translations[i] = 0.0;
        
        if (globals->trans_flags.estimate_center)
          globals->trans_info.center[i]       = c1[i+1];
        else 
          {
            if (!center_forced)
              globals->trans_info.center[i]       = 0.0;
          }
        if (globals->trans_flags.estimate_scale)
          globals->trans_info.scales[i]       = sc[i+1]; 
        else
          globals->trans_info.scales[i]       = 1.0;
        
        if (globals->trans_flags.estimate_quats)
          globals->trans_info.quaternions[i]   = quats[i];
        else
          globals->trans_info.quaternions[i]   = 0.0;
        
      }

      if (globals->trans_flags.estimate_quats)
        quats3=sqrt(1-SQR(globals->trans_info.quaternions[0])-SQR(globals->trans_info.quaternions[1])-SQR(globals->trans_info.quaternions[2]));
    
      VIO_FREE2D(trans);
      VIO_FREE2D(rots);
      FREE(quats);
      FREE(ang);
      FREE(c1);
      FREE(c2);
      FREE(sc);
  
   
	  if (globals->flags.debug) {
    
	    print ( "Transform center   = %8.3f %8.3f %8.3f\n", 
	            globals->trans_info.center[0],
	            globals->trans_info.center[1],
	            globals->trans_info.center[2] );
	    if(globals->trans_info.rotation_type == TRANS_ROT)
	      print ( "Transform rots    = %8.3f %8.3f %8.3f\n", 
	              globals->trans_info.rotations[0],
	              globals->trans_info.rotations[1],
	              globals->trans_info.rotations[2] );
	    else
	      print ( "Transform quaternion   = %8.3f %8.3f %8.3f %8.3f \n\n", 
	              globals->trans_info.quaternions[0],
	              globals->trans_info.quaternions[1],
	              globals->trans_info.quaternions[2],
	              quats3 );
    
	    print ( "Transform trans    = %8.3f %8.3f %8.3f\n", 
	            globals->trans_info.translations[0],
	            globals->trans_info.translations[1],
	            globals->trans_info.translations[2] );
	    print ( "Transform scale    = %8.3f %8.3f %8.3f\n\n", 
	            globals->trans_info.scales[0],
	            globals->trans_info.scales[1],
	            globals->trans_info.scales[2] );
	  }
	
  }
  

  else 
    { /*  we have an input matrix, we now have to extract the proper parameters from it */
      if (globals->flags.debug) print ("  using input transformation to get initial parameters:\n");
      	
      if (get_transform_type(globals->trans_info.transformation) == CONCATENATED_TRANSFORM) {
        lt = get_linear_transform_ptr(get_nth_general_transform(globals->trans_info.transformation,0));
      }
      else
        lt = get_linear_transform_ptr(globals->trans_info.transformation);

      
      /* get cog of data1 before extracting parameters
         from matrix, if estimate requested on command line */
      
      
      /* if centroid not forced on command
         line, then set it to 0,0,0 */
      
      if  (globals->trans_info.center[0] == -DBL_MAX &&
           globals->trans_info.center[1] == -DBL_MAX &&
           globals->trans_info.center[2] == -DBL_MAX) 
        {
          for(i=0; i<3; i++) globals->trans_info.center[i] = 0.0;
          
          if (globals->flags.debug) 
            {
              print ("   Center of rot/scale not forced, will be set to : %f %f %f\n",
                     globals->trans_info.center[0],
                     globals->trans_info.center[1],
                     globals->trans_info.center[2]);
            }
          
        }

      if (globals->trans_flags.estimate_center) 
        {  


          VIO_ALLOC2D(cov1 ,4,4);
          ALLOC(c1   ,4);
          
          if (! vol_cog(d1, m1,  c1, globals->step ) ) 
            {
              print_error_and_line_num("%s", __FILE__, __LINE__,"Cannot calculate the COG of volume 1\n." );
              return(FALSE);
            }

          for(i=0; i<=2; i++)
            globals->trans_info.center[i] = c1[i+1];
          
          if (globals->flags.debug) 
            print ("   User-requested COG estimate: %f %f %f\n",c1[1],c1[2],c1[3]);
          
          VIO_FREE2D(cov1);
          FREE(c1);
        }
                                               /* get the parameters from the input matrix: */
      if(globals->trans_info.rotation_type == TRANS_ROT)
        if (!extract2_parameters_from_matrix(lt,
                                             globals->trans_info.center,
                                             globals->trans_info.translations,
                                             globals->trans_info.scales,
                                             globals->trans_info.shears,
                                             globals->trans_info.rotations)) {
          return(FALSE);  
        }
      if(globals->trans_info.rotation_type == TRANS_QUAT)
        if (!extract2_parameters_from_matrix_quater(lt,
                                                    globals->trans_info.center,
                                                    globals->trans_info.translations,
                                                    globals->trans_info.scales,
                                                    globals->trans_info.shears,
                                                    globals->trans_info.quaternions)) {
          return(FALSE);  
        }
      
                                        /* do we need to replace anything?
                                           (note that center was possibly replaced
                                           just abve, and we dont have to do all the 
                                           PAT stuff if nothing else is needed)    */
      if ((globals->trans_flags.estimate_scale  ||
           globals->trans_flags.estimate_rots   ||
           globals->trans_flags.estimate_quats  ||
           globals->trans_flags.estimate_trans)) {
        
        if (globals->flags.debug) {
          print ("  using PAT to get overiding parameters:\n");
          print ("  will try to get:");
          if (globals->trans_flags.estimate_center) print (" [center already done]");
          if (globals->trans_flags.estimate_scale)  print (" [scale]");
          if (globals->trans_flags.estimate_rots)   print (" [rots]");
          if (globals->trans_flags.estimate_quats)   print (" [quats]");
          if (globals->trans_flags.estimate_trans)  print (" [trans]");
          print ("\n");
        }
        
        VIO_ALLOC2D(trans,5,5);
        VIO_ALLOC2D(rots,5,5);
        ALLOC(quats ,4);
        ALLOC(ang ,4);
        ALLOC(c1 ,4);
        ALLOC(c2 ,4);
        ALLOC(sc ,4);
        
        for(i=0; i<3; i++) c1[i+1] = globals->trans_info.center[i];
        
        if(globals->trans_info.rotation_type == TRANS_ROT)
          if (!init_transformation(d1,d2,m1,m2, globals->step, globals->flags.verbose,
                                   trans,rots,ang,c1,c2,sc, 
                                   center_forced,&(globals->trans_flags)))
            return(FALSE);      
        if(globals->trans_info.rotation_type == TRANS_QUAT)
          if (!init_transformation_quater(d1,d2,m1,m2, globals->step, globals->flags.verbose,
                                          trans,ang,quats,c1,c2,sc, 
                                          center_forced,&(globals->trans_flags)))
            return(FALSE);      
      
        for(i=0; i<3; i++) 
          {
            if (globals->trans_flags.estimate_rots)
              globals->trans_info.rotations[i]    = ang[i+1];
            else
              globals->trans_info.rotations[i]    = 0.0;
            
            if (globals->trans_flags.estimate_trans)
              globals->trans_info.translations[i] = trans[i+1][4];
            else
              globals->trans_info.translations[i] = 0.0;
            
            if (globals->trans_flags.estimate_scale)
              globals->trans_info.scales[i]       = sc[i+1]; 
            else
              globals->trans_info.scales[i]       = 1.0;
            if (globals->trans_flags.estimate_quats)
              globals->trans_info.quaternions[i]   = quats[i]; 
            else
              globals->trans_info.quaternions[i]     = 0.0;
          }
   
      
        VIO_FREE2D(trans);
        VIO_FREE2D(rots);
        FREE(quats);
        FREE(ang);
        FREE(c1);
        FREE(c2);
        FREE(sc);
      }
      
    }
  

  if (get_transform_type(globals->trans_info.transformation) == CONCATENATED_TRANSFORM) {
    lt = get_linear_transform_ptr(get_nth_general_transform(globals->trans_info.transformation,0));
  }
  else {
    lt = get_linear_transform_ptr(globals->trans_info.transformation);
}

  if( globals->trans_info.rotation_type == TRANS_ROT)
    build_transformation_matrix(lt,
                                globals->trans_info.center,
                                globals->trans_info.translations,
                                globals->trans_info.scales,
                                globals->trans_info.shears,
                                globals->trans_info.rotations);

  if( globals->trans_info.rotation_type == TRANS_QUAT)
    {
      globals->trans_info.quaternions[3]=sqrt(1-SQR(globals->trans_info.quaternions[0])-SQR(globals->trans_info.quaternions[1])-SQR(globals->trans_info.quaternions[2]));
      build_transformation_matrix_quater(lt,
                                         globals->trans_info.center,
                                         globals->trans_info.translations,
                                         globals->trans_info.scales,
                                         globals->trans_info.shears,
                                         globals->trans_info.quaternions);
    }

  return(TRUE);
}






