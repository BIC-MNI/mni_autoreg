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
@MODIFIED   : Revision 1.5  1995-02-22 08:56:06  louis
@MODIFIED   : Montreal Neurological Institute version.
@MODIFIED   : compiled and working on SGI.  this is before any changes for SPARC/
@MODIFIED   : Solaris.
@MODIFIED   :
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
#include <recipes.h>

#include "constants.h"
#include "matrix_basics.h"
#include "cov_to_praxes.h"
#include "make_rots.h"
#include "point_vector.h"

#define INTERPOLATE_TRUE_VALUE(volume, coord, result) \
   trilinear_interpolant(volume, coord, result)

static char *default_dim_names[N_DIMENSIONS] =
   { MIzspace, MIyspace, MIxspace };



char *prog_name;


BOOLEAN vol_to_cov(Volume d1, Volume m1, float *centroid, float **covar, double *step)
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

  Real
    tx,ty,tz,
    voxel_value;
  int
    i,r,c,s,
    limits[VOL_NDIMS];

  float
    t,
    sxx,syy,szz,
    sxy,syz,sxz,
    sx,sy,sz,si; 
  Real
    thickness[3];
  int
    sizes[3];

  Real
    true_value,
    position[3];

  get_volume_separations(d1, thickness);
  get_volume_sizes(d1, sizes);
  
				/* build sampling lattice info */
  for_less( i, 0, 3) {	
    step[i] *= thickness[i] / ABS( thickness[i]);
  }

  fill_Vector( col_step,   step[COL_IND], 0.0,     0.0 );
  fill_Vector( row_step,   0.0,     step[ROW_IND], 0.0 );
  fill_Vector( slice_step, 0.0,     0.0,     step[SLICE_IND] );

  convert_3D_voxel_to_world(d1, 0.0, 0.0, 0.0, &tx, &ty, &tz); 

  fill_Point( starting_origin, tx, ty, tz);

  for_less( i, 0, 3) {		/* for each dim, get # of steps in that direction,
				   and set starting offset */
    t = sizes[i] * thickness[i] / step[i];
    limits[i] = (int)( ABS( t ) );
    
    Point_coord( starting_offset, (i) ) = 
      ( (sizes[i]-1)*thickness[i] - (limits[i] * step[i] ) ) / 2.0;
  }
  
  ADD_POINTS( starting_position, starting_origin, starting_offset ); /*  */

				/* calculate centroids first */

  sx = 0.0;
  sy = 0.0;
  sz = 0.0;
  si = 0.0;

  for_inclusive(s,0,limits[SLICE_IND]) {

    SCALE_VECTOR( vector_step, slice_step, s);
    ADD_POINT_VECTOR( slice, starting_position, vector_step );

    for_inclusive(r,0,limits[ROW_IND]) {

      SCALE_VECTOR( vector_step, row_step, r);
      ADD_POINT_VECTOR( row, slice, vector_step );

      SCALE_POINT( col, row, 1.0); /* init first col position */
      for_inclusive(c,0,limits[COL_IND]) {

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

    for_inclusive(s,0,limits[Z]) {
      
      SCALE_VECTOR( vector_step, slice_step, s);
      ADD_POINT_VECTOR( slice, starting_position, vector_step );
      
      for_inclusive(r,0,limits[Y]) {
	
	SCALE_VECTOR( vector_step, row_step, r);
	ADD_POINT_VECTOR( row, slice, vector_step );
	
	SCALE_POINT( col, row, 1.0); /* init first col position */
	for_inclusive(c,0,limits[X]) {
	  
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



BOOLEAN get_cog(char *file, double *c1)
{
  Volume vol;
  float **cov, *cog;
  double step[3];
  Real x,y,z,r,s,c;

  input_volume(file,3,default_dim_names /*(char **)NULL*/, NC_UNSPECIFIED, FALSE, 0.0,0.0,
	       TRUE, &vol, (minc_input_options *)NULL);

  step[0] = 4.0;
  step[1] = 4.0;
  step[2] = 4.0;

  cov = matrix(1,3,1,3); 
  cog = vector(1,3);

  if ( vol_to_cov(vol, NULL, cog, cov, step ) ) {
    c1[0] = cog[1];
    c1[1] = cog[2];
    c1[2] = cog[3];
    return(TRUE);
  }
  else
    return(FALSE);


}




/* Main program */

int main(int argc, char *argv[])
{
   General_transform transform, new_transform;
   Transform
     *lt;
   Real scales[3], trans[3], rots[3], skews[3], center[3];
   int i;

   prog_name = argv[0];

   /* Check arguments */
   if (argc != 3 && argc != 4) {
      (void) fprintf(stderr, "Usage: %s <input.xfm> <result.xfm> [<file.mnc>]\n",
                     argv[0]);
      exit(EXIT_FAILURE);
   }

   /* Read in file to check scale */
   if (input_transform_file(argv[1], &transform) != OK) {
      (void) fprintf(stderr, "%s: Error reading transform file %s\n",
                     argv[0], argv[1]);
      exit(EXIT_FAILURE);
   }

   for_less(i,0,3)
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
   if (output_transform_file(argv[2], NULL, &new_transform) != OK) {
      (void) fprintf(stderr, "%s: Error writing transform file %s\n",
                     argv[0], argv[2]);
      exit(EXIT_FAILURE);
   }

   exit(EXIT_SUCCESS);
}
