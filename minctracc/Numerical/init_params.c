#include <def_mni.h>
#include <recipes.h>
#include "minctracc.h"
#include "def_geometry.h"

#include "matrix_basics.h"
#include "cov_to_praxes.h"
#include "make_rots.h"

extern Arg_Data main_args;

extern void 
  print_error(char *s, char * d1, int d2, int d3, int d4, int d5, int d6, int d7);

/* ----------------------------- MNI Header -----------------------------------
@NAME       : init_params.c
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: collection of routines that will calculate the parameters necessary
              from an input transformation matrix for optimization when 
	      mapping worls coordinates  of volume 1 into world coords in volume 2.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Thu May 27 16:50:50 EST 1993
                  
@MODIFIED   : 
---------------------------------------------------------------------------- */

/* ----------------------------- MNI Header -----------------------------------
@NAME       : vol_to_cov - get covariance and cog of volume.
@INPUT      : d1:
                one volume of data (already in memory).
	      kernel_size:
                size of resampling kernel
	        
@OUTPUT     : centroid - vector giving centroid of points. This vector
                         must be defined by the calling routine.
	      covar    - covariance matrix (in numerical recipes form).
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
public Boolean vol_to_cov(Volume d1, float *centroid, float **covar, Arg_Data *globals)
{

  Vector
    vector_step,
    slice_step,
    row_step,
    col_step;

  Point 
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

  get_volume_separations(d1, thickness);
  get_volume_sizes(d1, sizes);
  
				/* build sampling lattice info */
  for_less( i, 0, 3) {	
    globals->step[X] *= thickness[i] / ABS( thickness[i]);
  }

  fill_Vector( row_step,   globals->step[X], 0.0,           0.0 );
  fill_Vector( col_step,   0.0,           globals->step[Y], 0.0 );
  fill_Vector( slice_step, 0.0,           0.0,           globals->step[Z] );

  convert_voxel_to_world(d1, 0.0, 0.0, 0.0, &tx, &ty, &tz);

  fill_Point( starting_origin, tx, ty, tz);

  for_less( i, 0, 3) {		/* for each dim, get # of steps in that direction,
				   and set starting offset */
    t = sizes[i] * thickness[i] / globals->step[i];
    limits[i] = (int)( ABS( t ) );
    
    Point_coord( starting_offset, i ) = 
      ( (sizes[i]-1)*thickness[i] - (limits[i] * globals->step[i] ) ) / 2.0;
  }
  
  ADD_POINTS( starting_position, starting_origin, starting_offset ); /*  */
  


				/* calculate centroids first */

  sx = 0.0;
  sy = 0.0;
  sz = 0.0;
  si = 0.0;

  for_inclusive(s,0,limits[Z]) {

    SCALE_VECTOR( vector_step, slice_step, s);
    ADD_POINT_VECTOR( slice, starting_position, vector_step );

    for_inclusive(r,0,limits[Y]) {

      SCALE_VECTOR( vector_step, row_step, r);
      ADD_POINT_VECTOR( row, slice, vector_step );

      SCALE_POINT( col, row, 1.0); /* init first col position */
      for_inclusive(c,0,limits[X]) {

	convert_world_to_voxel(d1, Point_x(col), Point_y(col), Point_z(col), &tx, &ty, &tz);

	fill_Point( voxel, tx, ty, tz );

/*	trilinear_interpolant( d1, &voxel, &voxel_value );*/

	INTERPOLATE_VOXEL_VALUE( d1, &voxel, &voxel_value ); 
	
	sx +=  Point_x(col) * voxel_value;
	sy +=  Point_y(col) * voxel_value;
	sz +=  Point_z(col) * voxel_value;

	si += voxel_value;

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
	  
	  convert_world_to_voxel(d1, Point_x(col), Point_y(col), Point_z(col), &tx, &ty, &tz);

	  fill_Point( voxel, tx, ty, tz );
	  
	  trilinear_interpolant( d1, &voxel, &voxel_value );
	  
/*	  INTERPOLATE_VOXEL_VALUE( d1, &voxel, &voxel_value );*/
	  
	  sxx += (Point_x( col )-centroid[1]) * (Point_x( col )-centroid[1]) * voxel_value;
	  syy += (Point_y( col )-centroid[2]) * (Point_y( col )-centroid[2]) * voxel_value;
	  szz += (Point_z( col )-centroid[3]) * (Point_z( col )-centroid[3]) * voxel_value;
	  sxy += (Point_x( col )-centroid[1]) * (Point_y( col )-centroid[2]) * voxel_value;
	  syz += (Point_y( col )-centroid[2]) * (Point_z( col )-centroid[3]) * voxel_value;
	  sxz += (Point_x( col )-centroid[1]) * (Point_z( col )-centroid[3]) * voxel_value;

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



/* ----------------------------- MNI Header -----------------------------------
@NAME       : init_transformation - get trans parameters using principal axis
                 transformation technique.  the transformation points in 
		 d1 -> d2.
@INPUT      : d1,d2:
                two volumes of data (already in memory).
	      m1,m2:
                two mask volumes for data (already in memory).
	      main_args: globval argument structure
	      kx,ky,kz:
                size of resampling kernel
	        
@OUTPUT     : c1,c2    - vector giving volume centroids. This vector
                         must be defined by the calling routine.
	      trans    - translation matrix (in numerical recipes form).
                         must be defined by the calling routine.
	      rots     - rotation matrix (in numerical recipes form).
                         must be defined by the calling routine.
	      ang      - vector of rotation angles, in radians (in numerical recipes form).
                         must be defined by the calling routine.
	      c1       - vector for centroid of d1 (in numerical recipes form).
                         must be defined by the calling routine.
	      c2       - vector for centroid of d2 (in numerical recipes form).
                         must be defined by the calling routine.
@RETURNS    : TRUE if ok, FALSE if error.
@DESCRIPTION: this routine uses the principal axis method to determine the 
              initial guess for registration transformation
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb 5, 1992 lc
@MODIFIED   : Thu May 27 16:50:50 EST 1993 lc
                 rewrite for minc files and david's library
---------------------------------------------------------------------------- */
private  Boolean init_transformation(
				  Volume d1, /* data for volume1 */
				  Volume d2, /* data for volume2 */
				  Volume m1, /* mask for volume1 */
				  Volume m2, /* mask for volume2 */
				  Arg_Data *globals,
				  float kx,          /* size of current resampling kernel */
				  float ky,
				  float kz,
				  float **trans,     /* translation matrix to go from d1 to d2 */
				  float **rots,      /* rotation matrix to go from d1 to d2    */
				  float *ang,        /* rotation angles to go from d1 to d2    */
				  float *c1,         /* centroid of masked d1 */
				  float *c2,         /* centroid of masked d1 */
				  float *scale)      /* scaling from d1 to d2 */
{
  float
    tx,ty,tz,
    *angles,			/* rotation angles - rx,ry and rz */
    **cov1,**cov2,		/* covariance matrix */
    **prin_axes1, **prin_axes2, /* principal axis */
    **R1,**R2,			/* rotation matrix (normalized prin_axes) */
    **Rinv,**R;
  
  float
    kern,norm;
  
  int
    ndim,i,j;
  
  nr_identf(trans,1,4,1,4);	/* start with identity                       */
  nr_identf(rots, 1,4,1,4);
  
  cov1       = matrix(1,3,1,3); 
  cov2       = matrix(1,3,1,3); 
  prin_axes1 = matrix(1,3,1,3); 
  prin_axes2 = matrix(1,3,1,3); 
  R1         = matrix(1,3,1,3); 
  R2         = matrix(1,3,1,3); 
  R          = matrix(1,3,1,3); 
  Rinv       = matrix(1,3,1,3); 
  angles     = vector(1,3);
  
  
  kern = MIN( kx, ky);
  kern = MIN( kern, kz);
  
  if (! vol_to_cov(d1, c1, cov1, globals ) ) {
    print_error("Cannot calculate the COG of volume 1\n.", __FILE__, __LINE__, 0,0,0,0,0 );
    return(FALSE);
  }
  if (! vol_to_cov(d2, c2, cov2, globals ) ) {
    print_error("Cannot calculate the COG of volume 2\n.", __FILE__, __LINE__, 0,0,0,0,0 );
    return(FALSE);
  }

  tx = c2[1] - c1[1];    /* translations to map vol1 into vol2                  */
  ty = c2[2] - c1[2];
  tz = c2[3] - c1[3];
  
  cov_to_praxes(3, cov1, prin_axes1);   
  cov_to_praxes(3, cov2, prin_axes2);
  
  nr_copyf(prin_axes1,1,3,1,3,R1);
  nr_copyf(prin_axes2,1,3,1,3,R2);
  
  
  if (globals->flags.verbose > 1) {
    printf ("cov1:                              cov2:\n");
    for (i=1; i<=3; i++) {
      for (j=1; j<=3; j++)
	printf ("%8.3f ", cov1[i][j]);
    printf ("|");
      for (j=1; j<=3; j++)
	printf ("%8.3f ", cov2[i][j]);
      printf ("\n\n");
    }
  }

  for (i=1; i<=3; ++i) {
    norm = fsqrt(prin_axes1[i][1]*prin_axes1[i][1] + 
		 prin_axes1[i][2]*prin_axes1[i][2] + 
		 prin_axes1[i][3]*prin_axes1[i][3]);
    for (j=1; j<=3; ++j)
      R1[i][j] /= norm;
    
    norm = fsqrt(prin_axes2[i][1]*prin_axes2[i][1] + 
		 prin_axes2[i][2]*prin_axes2[i][2] + 
		 prin_axes2[i][3]*prin_axes2[i][3]);
    for (j=1; j<=3; ++j)
      R2[i][j] /= norm;
  }
  
  if (globals->flags.verbose > 1) {
    printf ("prin_axes1:                      princ_axes2:\n");
    for (i=1; i<=3; i++) {
      for (j=1; j<=3; j++)
	printf ("%8.3f ", prin_axes1[i][j]);
      printf ("|");
      for (j=1; j<=3; j++)
	printf ("%8.3f ", prin_axes2[i][j]);
      printf ("\n\n");
    }
  }
  
  invertmatrix(3,R1,Rinv);

  nr_multf(Rinv,1,3,1,3, R2,1,3,1,3, R);
  
  
  if (globals->flags.verbose > 1) {
    printf ("r1:                   r2:                  r:\n");
    for (i=1; i<=3; i++) {
      for (j=1; j<=3; j++)
	printf ("%8.3f ", R1[i][j]);
      printf ("|");
      for (j=1; j<=3; j++)
	printf ("%8.3f ", R2[i][j]);
      printf ("|");
      for (j=1; j<=3; j++)
	printf ("%8.3f ", R[i][j]);
      printf ("\n");
    }
  }


  transpose(3,3,R,R);		/* all of the princ axes stuff uses vec*mat */

  if (!rotmat_to_ang(R, angles)) {
    fprintf(stderr,"Could not extract angles from:\n");
    printmatrix(3,3,R);
    return(FALSE);
  }
  
  scale[1] = 1.0;
  scale[2] = 1.0;
  scale[3] = 1.0;

  ang[1] = angles[1];		/* rotation about X axis                         */
  ang[2] = angles[2];		/* rotation about Y axis                         */
  ang[3] = angles[3];		/* rotation about Z axis                         */
  
  for (i=1; i<=3; ++i)		/* set rotations in matrix                       */
    for (j=1; j<=3; ++j) {
      rots[i][j] = R[i][j];
    }
  
  trans[1][4] += tx;		/* set translations in translation matrix        */
  trans[2][4] += ty;
  trans[3][4] += tz;
  
  ndim = 3;
  
  if (globals->flags.verbose > 1) {
    (void) printf("\nFor volume 1 :");
    (void) printf("\nCentroid :");
    for (i=1; i<=ndim; i++) (void) printf("  %f",c1[i]);
    (void) printf("\n\n");
    (void) printf("Principal axes\n");
    for (i=1; i<=ndim; i++) {
      (void) printf("Vector %d :",i);
      for (j=1; j<=ndim; j++) {
	(void) printf("  %f",prin_axes1[i][j]);
      }
      (void) printf("\n");
    }
  
    (void) printf("\n");
    
    (void) printf("\nFor volume 2 :");
    (void) printf("\nCentroid :");
    for (i=1; i<=ndim; i++) (void) printf("  %f",c2[i]);
    (void) printf("\n\n");
    (void) printf("Principal axes\n");
    for (i=1; i<=ndim; i++) {
      (void) printf("Vector %d :",i);
      for (j=1; j<=ndim; j++) {
	(void) printf("  %f",prin_axes2[i][j]);
      }
      (void) printf("\n");
    }
    
    (void) printf ("rotation angles are: %f %f %f\n",angles[1]*180/PI,angles[2]*180/PI,
		   angles[3]*180/PI);
  }

  free_matrix(cov1,      1,3,1,3); 
  free_matrix(cov2,      1,3,1,3); 
  free_matrix(prin_axes1,1,3,1,3); 
  free_matrix(prin_axes2,1,3,1,3); 
  free_matrix(R1,        1,3,1,3); 
  free_matrix(R2,        1,3,1,3); 
  free_matrix(R,         1,3,1,3); 
  free_matrix(Rinv,      1,3,1,3); 
  free_vector(c1,        1,3);
  free_vector(c2,        1,3);
  free_vector(angles,    1,3);
  
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

public Boolean init_params(Volume d1,
		 Volume d2,
		 Volume m1,
		 Volume m2, 
		 Arg_Data *globals)
{  
  float 
    **trans,			/* principal componant variables */
    **rots,
    **cov1,
    *c1,
    *c2,
    *ang,
    *sc;
    
  int
    i;
    
  Linear_Transformation 
    *lt;

				/* if no transformation specified on the command line,
				   of if PAT tgransformation selected,
				     do principal axes transformation to init trans. */

  if (globals->trans_info.use_default == TRUE || 
      globals->trans_info.transform_type==TRANS_PAT) {

				/* if none of these flags set on the command line,
				   or if PAT transformation requested, 
				   then assume them all to be true, otherwise leave
				   them as set on the command line. */
    
    if ((globals->trans_info.transform_type==TRANS_PAT) || 
	!(globals->trans_flags.estimate_center ||
	  globals->trans_flags.estimate_scale  ||
	  globals->trans_flags.estimate_rots   ||
	  globals->trans_flags.estimate_trans)) {
      
      globals->trans_flags.estimate_center = TRUE;
      globals->trans_flags.estimate_scale = TRUE;
      globals->trans_flags.estimate_rots = TRUE;
      globals->trans_flags.estimate_trans = TRUE;
    }
    
    trans = matrix(1,4,1,4);
    rots = matrix(1,4,1,4);
    ang = vector(1,3);
    c1 = vector(1,3);
    c2 = vector(1,3);
    sc = vector(1,3);
    
    if (!init_transformation(d1,d2,m1,m2, globals,
			     globals->step[0],globals->step[1],globals->step[2],
			     trans,rots,ang,c1,c2,sc))
      return(FALSE);
    
    for_less( i, 0, 3 ) {
      if (globals->trans_flags.estimate_rots)
	globals->trans_info.rotations[i]    = ang[i+1]*RAD_TO_DEG;
      else
	globals->trans_info.rotations[i]    = 0.0;
      
      if (globals->trans_flags.estimate_trans)
	globals->trans_info.translations[i] = trans[i+1][4];
      else
	globals->trans_info.translations[i] = 0.0;
      
      if (globals->trans_flags.estimate_center)
	globals->trans_info.center[i]       = c1[i+1];
      else
	globals->trans_info.center[i]       = 0.0;
      
      if (globals->trans_flags.estimate_scale)
	globals->trans_info.scales[i]       = sc[i+1]; 
      else
	globals->trans_info.scales[i]       = 1.0;
    }
    
    free_matrix(trans,1,4,1,4);
    free_matrix(rots,1,4,1,4);
    free_vector(ang,1,3);
    free_vector(c1,1,3);
    free_vector(c2,1,3);
    free_vector(sc,1,3);
    
  }
  else { /*  we have an input matrix, we now have to extract the proper parameters from it */

    lt = (Linear_Transformation *) globals->trans_info.transformation.trans_data;

				/* get cog of data1 before extracting parameters
				   from matrix, if estimate requested on command line */
    if (globals->trans_flags.estimate_center) {  

      cov1 = matrix(1,3,1,3); 
      c1   = vector(1,3);
      
      if (! vol_to_cov(d1, c1, cov1, globals ) ) {
	print_error("Cannot calculate the COG of volume 1\n.", 
		    __FILE__, __LINE__, 0,0,0,0,0 );
	return(FALSE);
      }
      for_inclusive( i, 0, 2 )
	globals->trans_info.center[i] = c1[i+1];

      free_matrix(cov1,1,3,1,3);
      free_vector(c1,1,3);
    }
				/* get the parameters from the input matrix: */
    extract_parameters_from_matrix(lt->mat, 
				   globals->trans_info.center,
				   globals->trans_info.translations,
				   globals->trans_info.scales,
				   globals->trans_info.rotations);

				/* do we need to replace anything?
				   (note that center was possibly replaced
				    just abve, and we dont have to do all the 
				    PAT stuff if nothing else is needed)         */
    if ((globals->trans_flags.estimate_scale  ||
	 globals->trans_flags.estimate_rots   ||
	 globals->trans_flags.estimate_trans)) {
      
      trans = matrix(1,4,1,4);
      rots = matrix(1,4,1,4);
      ang = vector(1,3);
      c1 = vector(1,3);
      c2 = vector(1,3);
      sc = vector(1,3);
    
      init_transformation(d1,d2,m1,m2, globals,
			  globals->step[0],globals->step[1],globals->step[2],
			  trans,rots,ang,c1,c2,sc);      
      
      for_less( i, 0, 3 ) {
	if (globals->trans_flags.estimate_rots)
	  globals->trans_info.rotations[i]    = ang[i+1]*RAD_TO_DEG;
	else
	  globals->trans_info.rotations[i]    = 0.0;
	
	if (globals->trans_flags.estimate_trans)
	  globals->trans_info.translations[i] = trans[i+1][4];
	else
	  globals->trans_info.translations[i] = 0.0;
	
	if (globals->trans_flags.estimate_center)
	  globals->trans_info.center[i]       = c1[i+1];
	else
	  globals->trans_info.center[i]       = 0.0;
	
	if (globals->trans_flags.estimate_scale)
	  globals->trans_info.scales[i]       = sc[i+1]; 
	else
	  globals->trans_info.scales[i]       = 1.0;
      }
        
    free_matrix(trans,1,4,1,4);
    free_matrix(rots,1,4,1,4);
    free_vector(ang,1,3);
    free_vector(c1,1,3);
    free_vector(c2,1,3);
    free_vector(sc,1,3);
    }

  }

  lt = (Linear_Transformation *) globals->trans_info.transformation.trans_data;
  
  build_transformation_matrix(lt->mat, 
			      globals->trans_info.center,
			      globals->trans_info.translations,
			      globals->trans_info.scales,
			      globals->trans_info.rotations);
  

  return(TRUE);
}




