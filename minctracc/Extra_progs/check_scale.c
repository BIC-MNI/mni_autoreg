/* ----------------------------- MNI Header -----------------------------------
@NAME       : check_scale.c
@INPUT      : argc, argv - command line arguments
@OUTPUT     : (none)
@RETURNS    : status
@DESCRIPTION: Program to check the z-scale of an  MNI transform file
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Mon Nov 29 11:01:47 EST 1993 Louis
@MODIFIED   : $Log: check_scale.c,v $
@MODIFIED   : Revision 96.6  2006-11-30 09:07:31  rotor
@MODIFIED   :  * many more changes for clean minc 2.0 build
@MODIFIED   :
@MODIFIED   : Revision 96.5  2006/11/29 09:09:30  rotor
@MODIFIED   :  * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   : Revision 96.4  2005/07/20 20:45:45  rotor
@MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :     * Removed old VOLUME_IO cruft #defines
@MODIFIED   :     * Fixed up all Makefile.am's in subdirs
@MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   : Revision 96.3  2002/03/26 14:15:27  stever
@MODIFIED   : Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   : Revision 96.2  2000/03/15 08:42:36  stever
@MODIFIED   : Code cleanup: all functions prototyped (except ParseArgs.c), no useless declarations, etc
@MODIFIED   :
@MODIFIED   : Revision 96.1  1999/10/25 19:52:05  louis
@MODIFIED   : final checkin before switch to CVS
@MODIFIED   :
 * Revision 96.0  1996/08/21  18:21:36  louis
 * Release of MNI_AutoReg version 0.96
 *
 * Revision 9.6  1996/08/21  18:21:30  louis
 * Pre-release
 *
 * Revision 9.5  1996/08/12  14:15:18  louis
 * Never released version 0.95
 *
 * Revision 1.6  1996/08/12  14:15:08  louis
 * Pre-release
 *
 * Revision 1.5  1995/02/22  08:56:06  collins
 * Montreal Neurological Institute version.
 * compiled and working on SGI.  this is before any changes for SPARC/
 * Solaris.
 *
 * Revision 1.4  94/06/06  09:30:58  louis
 * *** empty log message ***
 * 
 * Revision 1.3  94/04/26  12:52:49  louis
 * updated with new versions of make_rots, extract2_parameters_from_matrix 
 * that include proper interpretation of skew.
 * 
 * Revision 1.2  94/04/06  11:46:45  louis
 * working linted version of linear + non-linear registration based on Lvv
 * operator working in 3D
 * 
 * Revision 1.1  94/02/21  16:31:52  louis
 * Initial revision
 * 
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
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="";
#endif

#include <volume_io.h>

#include "constants.h"
#include "matrix_basics.h"
#include "make_rots.h"
#include "point_vector.h"
#include "interpolation.h"

#define INTERPOLATE_TRUE_VALUE(volume, coord, result) \
   trilinear_interpolant(volume, coord, result)

static char *default_dim_names[VIO_N_DIMENSIONS] =
   { MIzspace, MIyspace, MIxspace };



char *prog_name;


VIO_BOOL vol_to_cov(VIO_Volume d1, VIO_Volume m1, float *centroid, float **covar, double *step)
{

  VectorR
    vector_step,
    slice_step,
    row_step,
    col_step;

  PointR 
    starting_offset,
    starting_origin,
    starting_position,
    slice,
    row,
    col,
    voxel;

  VIO_Real
    tx,ty,tz;
  int
    i,r,c,s,
    limits[VOL_NDIMS];

  float
    t,
    sxx,syy,szz,
    sxy,syz,sxz,
    sx,sy,sz,si; 
  VIO_Real
    thickness[3];
  int
    sizes[3];

  VIO_Real
    true_value;

  get_volume_separations(d1, thickness);
  get_volume_sizes(d1, sizes);
  
                                /* build sampling lattice info */
  for(i=0; i<3; i++) {        
    step[i] *= thickness[i] / fabs( thickness[i]);
  }

  fill_Vector( col_step,   step[COL_IND], 0.0,     0.0 );
  fill_Vector( row_step,   0.0,     step[ROW_IND], 0.0 );
  fill_Vector( slice_step, 0.0,     0.0,     step[SLICE_IND] );

  convert_3D_voxel_to_world(d1, 0.0, 0.0, 0.0, &tx, &ty, &tz); 

  fill_Point( starting_origin, tx, ty, tz);

  for(i=0; i<3; i++) {                /* for each dim, get # of steps in that direction,
                                   and set starting offset */
    t = sizes[i] * thickness[i] / step[i];
    limits[i] = abs( t );
    
    Point_coord( starting_offset, (i) ) = 
      ( (sizes[i]-1)*thickness[i] - (limits[i] * step[i] ) ) / 2.0;
  }
  
  ADD_POINTS( starting_position, starting_origin, starting_offset ); /*  */

                                /* calculate centroids first */

  sx = 0.0;
  sy = 0.0;
  sz = 0.0;
  si = 0.0;

  for(s=0; s<=limits[SLICE_IND]; s++) {

    SCALE_VECTOR( vector_step, slice_step, s);
    ADD_POINT_VECTOR( slice, starting_position, vector_step );

    for(r=0; r<=limits[ROW_IND]; r++) {

      SCALE_VECTOR( vector_step, row_step, r);
      ADD_POINT_VECTOR( row, slice, vector_step );

      SCALE_POINT( col, row, 1.0); /* init first col position */
      for(c=0; c<=limits[COL_IND]; c++) {

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
        
        ADD_POINT_VECTOR( col, col, col_step );
        
      }
    }
  }

  if (si!=0.0) {
    centroid[1] = sx/ si;
    centroid[2] = sy/ si;
    centroid[3] = sz/ si;
    
    sxx = syy = szz = 0.0;
    sxy = syz = sxz = 0.0;
    
                                /* now calculate variances and co-variances */

    for(s=0; s<=limits[VIO_Z]; s++) {
      
      SCALE_VECTOR( vector_step, slice_step, s);
      ADD_POINT_VECTOR( slice, starting_position, vector_step );
      
      for(r=0; r<=limits[VIO_Y]; r++) {
        
        SCALE_VECTOR( vector_step, row_step, r);
        ADD_POINT_VECTOR( row, slice, vector_step );
        
        SCALE_POINT( col, row, 1.0); /* init first col position */
        for(c=0; c<=limits[VIO_X]; c++) {
          
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

            }
            /* else requested voxel is just outside volume., so ignore it */

          } 
          
          ADD_POINT_VECTOR( col, col, col_step );
        
        }
      }
    }
    
    covar[1][1] = sxx/si; covar[1][2] = sxy/si; covar[1][3] = sxz/si;
    covar[2][1] = sxy/si; covar[2][2] = syy/si; covar[2][3] = syz/si;
    covar[3][1] = sxz/si; covar[3][2] = syz/si; covar[3][3] = szz/si;
    
    return(TRUE);
    
  }
  else {
    return(FALSE);
  }
}



VIO_BOOL get_cog(char *file, double *c1)
{
  VIO_Volume vol;
  float **cov, *cog;
  double step[3];

  input_volume(file,3,default_dim_names /*(char **)NULL*/, NC_UNSPECIFIED, FALSE, 0.0,0.0,
               TRUE, &vol, (minc_input_options *)NULL);

  step[0] = 4.0;
  step[1] = 4.0;
  step[2] = 4.0;

  VIO_ALLOC2D(cov,4,4);
  ALLOC(cog,4);

  if ( vol_to_cov(vol, NULL, cog, cov, step ) ) {
    c1[0] = cog[1];
    c1[1] = cog[2];
    c1[2] = cog[3];
    return(TRUE);
  }
  else
    return(FALSE);

  VIO_FREE2D(cov);
  FREE(cog);

}




/* Main program */

int main(int argc, char *argv[])
{
   VIO_General_transform transform, new_transform;
   VIO_Transform
     *lt;
   VIO_Real scales[3], trans[3], rots[3], skews[3], center[3];
   int i;

   prog_name = argv[0];

   /* Check arguments */
   if (argc != 3 && argc != 4) {
      (void) fprintf(stderr, "Usage: %s <input.xfm> <result.xfm> [<file.mnc>]\n",
                     argv[0]);
      exit(EXIT_FAILURE);
   }

   /* Read in file to check scale */
   if (input_transform_file(argv[1], &transform) != VIO_OK) {
      (void) fprintf(stderr, "%s: Error reading transform file %s\n",
                     argv[0], argv[1]);
      exit(EXIT_FAILURE);
   }

   for(i=0; i<3; i++)
     center[i] = 0.0;
   if (argc>3) {
     print ("mnc = %s\n",argv[3]);
     if (! get_cog(argv[3], center) ) {
       print("Cannot calculate the COG of volume %s\n.", argv[3] );
       return(FALSE);
     }
   }


   if (get_transform_type(&transform) == CONCATENATED_TRANSFORM) {
     (void) fprintf(stderr, "Error: Cannot deal with concatenated transforms\n");
     exit(EXIT_FAILURE);
   }
   if (get_transform_type(&transform) == THIN_PLATE_SPLINE) {
     (void) fprintf(stderr, "Error: Cannot deal with non-linear transforms\n");
     exit(EXIT_FAILURE);
   }
   if (get_transform_type(&transform) == USER_TRANSFORM) {
     (void) fprintf(stderr, "Error: Cannot deal with user-defined transforms\n");
     exit(EXIT_FAILURE);
   }

   
   /* Extract parameters from transform */


   lt = get_linear_transform_ptr(&transform);

   if (!extract2_parameters_from_matrix(lt,
                                       center,
                                       trans,
                                       scales,
                                       skews,
                                       rots)) {
     (void) fprintf(stderr, "Error: Cannot extract parameters from matrix\n");
     exit(EXIT_FAILURE);
   }
        
   /* check scaling parameters */


   if (scales[2] > 1.15*(scales[0]+scales[1])/2.0) {
     scales[2] = (scales[0]+scales[1])/2.0;
   }
   
   /* rebuild transformation from parameters */

   build_transformation_matrix(lt,
                               center,
                               trans,
                               scales,
                               skews,
                               rots);

   create_linear_transform(&new_transform, lt);
   

   /* Write out the transform */
   if (output_transform_file(argv[2], NULL, &new_transform) != VIO_OK) {
      (void) fprintf(stderr, "%s: Error writing transform file %s\n",
                     argv[0], argv[2]);
      exit(EXIT_FAILURE);
   }

   exit(EXIT_SUCCESS);
}
