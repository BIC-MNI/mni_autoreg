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
@MODIFIED   :  Revision 1.10  1994-06-02 20:16:01  louis
@MODIFIED   :  made modifications to allow deformations to be calulated in 2D on slices.
@MODIFIED   :  changes had to be made in set_up_lattice, init_lattice when defining
@MODIFIED   :  the special case of a single slice....
@MODIFIED   :  Build_default_deformation_field also had to reflect these changes.
@MODIFIED   :  do_non-linear-optimization also had to check if one of dimensions had
@MODIFIED   :  a single element.
@MODIFIED   :  All these changes were made, and slightly tested.  Another type of
@MODIFIED   :  deformation strategy will be necessary (to replace the deformation
@MODIFIED   :  perpendicular to the surface, since it does not work well).
@MODIFIED   :
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
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Numerical/init_params.c,v 1.10 1994-06-02 20:16:01 louis Exp $";
#endif


#include <volume_io.h>
#include <recipes.h>

#include "constants.h"
#include "arg_data.h"

#include "matrix_basics.h"
#include "cov_to_praxes.h"
#include "make_rots.h"

extern Arg_Data main_args;

#include "local_macros.h"
#include <print_error.h>

public void set_up_lattice(Volume data, 
			   double *user_step, /* user requested spacing for lattice */
			   double *start,     /* world starting position of lattice */
			   int    *count,     /* number of steps in each direction */
			   double *step,      /* step size in each direction */
			   VectorR directions[]); /* array of vector directions for each index*/




/* ----------------------------- MNI Header -----------------------------------
@NAME       : vol_to_cov - get covariance and cog of volume.
@INPUT      : d1: one volume of data (already in memory).
	      m1: its corresponding mask volume
	      step: an 3 element array of step sizes in x,ya nd z directions
	        
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
BOOLEAN vol_to_cov(Volume d1, Volume m1, float *centroid, float **covar, double *step)
{


  VectorR
    vector_step;

  PointR 
    starting_position,
    slice,
    row,
    col,
    voxel;

  Real
    tx,ty,tz;

  int
    r,c,s;

  float
    sxx,syy,szz,
    sxy,syz,sxz,
    sx,sy,sz,si; 

  Real
    true_value;

  int
    i,count[3];
  double 
    start[3],
    local_step[3];
  VectorR
    directions[3];  

				/* build default sampling lattice info
				   on the data set (d1)               */
  set_up_lattice(d1, step,
		 start, count, local_step, directions);

  if (main_args.flags.debug) {
    print ("in vol to cov\n");
    print ("start = %8.2f %8.2f %8.2f \n",start[0],start[1],start[2]);
    print ("count = %8d %8d %8d \n",count[0],count[1],count[2]);
    print ("step  = %8.2f %8.2f %8.2f \n",local_step[0],local_step[1],local_step[2]);
    
    for_less(i,0,3)
      print ("direct= %8.2f %8.2f %8.2f \n",
	     Point_x(directions[i]),
	     Point_y(directions[i]),
	     Point_z(directions[i]));
  }


  fill_Point( starting_position, start[0], start[1], start[2]);
  
				/* calculate centroids first */

  sx = 0.0;
  sy = 0.0;
  sz = 0.0;
  si = 0.0;

  for_inclusive(s,0,count[SLICE_IND]) {

    SCALE_VECTOR( vector_step, directions[SLICE_IND], s);
    ADD_POINT_VECTOR( slice, starting_position, vector_step );

    for_inclusive(r,0,count[ROW_IND]) {

      SCALE_VECTOR( vector_step, directions[ROW_IND], r);
      ADD_POINT_VECTOR( row, slice, vector_step );

      SCALE_POINT( col, row, 1.0); /* init first col position */
      for_inclusive(c,0,count[COL_IND]) {

	
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
	
	ADD_POINT_VECTOR( col, col, directions[COL_IND] );
	
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

    for_inclusive(s,0,count[SLICE_IND]) {
      
      SCALE_VECTOR( vector_step, directions[SLICE_IND], s);
      ADD_POINT_VECTOR( slice, starting_position, vector_step );
      
      for_inclusive(r,0,count[ROW_IND]) {
	
	SCALE_VECTOR( vector_step, directions[ROW_IND], r);
	ADD_POINT_VECTOR( row, slice, vector_step );
	
	SCALE_POINT( col, row, 1.0); /* init first col position */
	for_inclusive(c,0,count[COL_IND]) {
	  
	  
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
	  
	  ADD_POINT_VECTOR( col, col, directions[COL_IND] );
	
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
	      step: an 3 element array of step sizes in x,y and z directions
	      verbose: =0 for quiet, >1 otherwise
	        
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
private  BOOLEAN init_transformation(
				     Volume d1, /* data for volume1 */
				     Volume d2, /* data for volume2 */
				     Volume m1, /* mask for volume1 */
				     Volume m2, /* mask for volume2 */
				     double *step, 
				     int    verbose,
				     float **trans,     /* translation matrix to go from d1 to d2 */
				     float **rots,      /* rotation matrix to go from d1 to d2    */
				     float *ang,        /* rotation angles to go from d1 to d2    */
				     float *c1,         /* centroid of masked d1 */
				     float *c2,         /* centroid of masked d1 */
				     float *scale,      /* scaling from d1 to d2 */
				     Transform_Flags *flags) /* flags for estimation */
{
  float
    tx,ty,tz,
    *angles,			/* rotation angles - rx,ry and rz */
    **cov1,**cov2,		/* covariance matrix */
    **prin_axes1, **prin_axes2, /* principal axis */
    **R1,**R2,			/* rotation matrix (normalized prin_axes) */
    **Rinv,**R;
  
  float
    norm;
  
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
  
  
  
  if (! vol_to_cov(d1, m1, c1, cov1, step ) ) {
    print_error("%s", __FILE__, __LINE__,"Cannot calculate the COG of volume 1\n." );
    return(FALSE);
  }

  if (verbose>0) print ("COG of v1: %f %f %f\n",c1[1],c1[2],c1[3]);

  if (flags->estimate_trans || flags->estimate_rots || flags->estimate_scale) {
    if (! vol_to_cov(d2, m2, c2, cov2, step ) ) {
      print_error("%s", __FILE__, __LINE__,"Cannot calculate the COG of volume 2\n." );
      return(FALSE);
    }
    if (verbose>0) print ("COG of v2: %f %f %f\n",c2[1],c2[2],c2[3]);
  }
  

  if (flags->estimate_trans) {
    tx = c2[1] - c1[1];    /* translations to map vol1 into vol2                  */
    ty = c2[2] - c1[2];
    tz = c2[3] - c1[3];
  }
  else {
    tx = ty = tz = 0.0;
  }

  trans[1][4] += tx;		/* set translations in translation matrix        */
  trans[2][4] += ty;
  trans[3][4] += tz;
  
  if (flags->estimate_rots || flags->estimate_scale) {
    cov_to_praxes(3, cov1, prin_axes1);   
    cov_to_praxes(3, cov2, prin_axes2);
    
    nr_copyf(prin_axes1,1,3,1,3,R1);
    nr_copyf(prin_axes2,1,3,1,3,R2);
    
    
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
    
    invertmatrix(3,R1,Rinv);
    
    nr_multf(Rinv,1,3,1,3, R2,1,3,1,3, R);
    
    
    if (verbose > 1) {
      print ("r1:                   r2:                  r:\n");
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
    
    
    transpose(3,3,R,R);		/* all of the princ axes stuff uses vec*mat */
    
    if (!rotmat_to_ang(R, angles)) {
      (void)fprintf(stderr,"Could not extract angles from:\n");
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
    
    ndim = 3;
    
    if (verbose > 1) {
      (void) print("\nFor volume 1 :");
      (void) print("\nCentroid :");
      for (i=1; i<=ndim; i++) (void) print("  %f",c1[i]);
      (void) print("\n\n");
      (void) print("Principal axes\n");
      for (i=1; i<=ndim; i++) {
	(void) print("Vector %d :",i);
	for (j=1; j<=ndim; j++) {
	  (void) print("  %f",prin_axes1[i][j]);
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
	(void) print("Vector %d :",i);
	for (j=1; j<=ndim; j++) {
	  (void) print("  %f",prin_axes2[i][j]);
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
    scale[1] = 1.0;
    scale[2] = 1.0;
    scale[3] = 1.0;
    
    ang[1] = 0.0;
    ang[2] = 0.0;
    ang[3] = 0.0;
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

public BOOLEAN init_params(Volume d1,
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
    
  Transform
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

    if (globals->filenames.measure_file!=NULL)

    if ((globals->trans_info.transform_type==TRANS_PAT) || 
	!(globals->filenames.measure_file!=NULL ||
	  globals->trans_flags.estimate_center ||
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
    
				/* get cog of data1 before extracting parameters
				   from matrix, if estimate requested on command line */



    if (!init_transformation(d1,d2,m1,m2, globals->step, globals->flags.verbose,
			     trans,rots,ang,c1,c2,sc, &(globals->trans_flags)))
      return(FALSE);
    
    for_less( i, 0, 3 ) {
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

    if (get_transform_type(globals->trans_info.transformation) == CONCATENATED_TRANSFORM) {
      lt = get_linear_transform_ptr(get_nth_general_transform(globals->trans_info.transformation,0));
    }
    else
      lt = get_linear_transform_ptr(globals->trans_info.transformation);

				/* get cog of data1 before extracting parameters
				   from matrix, if estimate requested on command line */
    if (globals->trans_flags.estimate_center) {  

      cov1 = matrix(1,3,1,3); 
      c1   = vector(1,3);
      
      if (! vol_to_cov(d1, m1,  c1, cov1, globals->step ) ) {
	print_error("%s", __FILE__, __LINE__,"Cannot calculate the COG of volume 1\n." );
	return(FALSE);
      }
      for_inclusive( i, 0, 2 )
	globals->trans_info.center[i] = c1[i+1];

      free_matrix(cov1,1,3,1,3);
      free_vector(c1,1,3);
    }
				/* get the parameters from the input matrix: */
    if (!extract2_parameters_from_matrix(lt,
					globals->trans_info.center,
					globals->trans_info.translations,
					globals->trans_info.scales,
					globals->trans_info.shears,
					globals->trans_info.rotations)) {
      return(FALSE);  
    }
	
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
    
      if (!init_transformation(d1,d2,m1,m2, globals->step, globals->flags.verbose,
 			       trans,rots,ang,c1,c2,sc, &(globals->trans_flags)))
	return(FALSE);      
      
      for_less( i, 0, 3 ) {
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

  if (get_transform_type(globals->trans_info.transformation) == CONCATENATED_TRANSFORM) {
    lt = get_linear_transform_ptr(get_nth_general_transform(globals->trans_info.transformation,0));
  }
  else
    lt = get_linear_transform_ptr(globals->trans_info.transformation);
  
    
  build_transformation_matrix(lt,
			      globals->trans_info.center,
			      globals->trans_info.translations,
			      globals->trans_info.scales,
			      globals->trans_info.shears,
			      globals->trans_info.rotations);
  

  return(TRUE);
}




