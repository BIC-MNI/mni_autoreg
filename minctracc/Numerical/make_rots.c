/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_rots.c
@DESCRIPTION: collection of routines to build transformation matrices
              from parameters and to extract parameters from matrices.
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
@CREATED    : Tue Jun  8 08:44:59 EST 1993 LC
@MODIFIED   : $Log: make_rots.c,v $
@MODIFIED   : Revision 1.7  1994-02-21 16:35:42  louis
@MODIFIED   : version before feb 22 changes
@MODIFIED   :
 * Revision 1.6  93/11/15  16:27:02  louis
 * working version, with new library, with RCS revision stuff,
 * before deformations included
 * 
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Numerical/make_rots.c,v 1.7 1994-02-21 16:35:42 louis Exp $";
#endif

#include <volume_io.h>
#include <recipes.h>

#include "matrix_basics.h"
#include "rotmat_to_ang.h"

#define  FILL_NR_COLVEC( vector, x, y, z ) \
            { \
		vector[1][1] = (x); \
		vector[2][1] = (y); \
		vector[3][1] = (z); \
		vector[4][1] = 1.0; \
            }

#define  MAG_NR_COLVEC( vector ) \
            ( fsqrt(vector[1][1]*vector[1][1] + \
		    vector[2][1]*vector[2][1] + \
		    vector[3][1]*vector[3][1])  \
	     )

#define  ADD_NR_COLVEC( vector, v1, v2 ) \
            { \
		vector[1][1] = v1[1][1] + v2[1][1]; \
		vector[2][1] = v1[2][1] + v2[2][1]; \
		vector[3][1] = v1[3][1] + v2[3][1]; \
		vector[4][1] = 1.0; \
            }

#define  SUB_NR_COLVEC( vector, v1, v2 ) \
            { \
		vector[1][1] = v1[1][1] - v2[1][1]; \
		vector[2][1] = v1[2][1] - v2[2][1]; \
		vector[3][1] = v1[3][1] - v2[3][1]; \
		vector[4][1] = 1.0; \
            }

/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_rots
@INPUT      : rot_x, rot_y, rot_z - three rotation angles, in radians.
@OUTPUT     : xmat, a numerical recipes matrix for homogeous transformations
@RETURNS    : nothing
@DESCRIPTION: to be applied by premultiplication, ie rot*vec = newvec
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Tue Jun  8 08:44:59 EST 1993 LC
@MODIFIED   : 
---------------------------------------------------------------------------- */

void   make_rots(float **xmat, float rot_x, float rot_y, float rot_z)
{
   float
      **TRX,
      **TRY,
      **TRZ,
      **T1,
      **T2;
   
   TRX  = matrix(1,4,1,4);
   TRY  = matrix(1,4,1,4);
   TRZ  = matrix(1,4,1,4);
   T1   = matrix(1,4,1,4);
   T2   = matrix(1,4,1,4);

   nr_rotxf(TRX, rot_x);             /* create the rotate X matrix */
   nr_rotyf(TRY, rot_y);             /* create the rotate Y matrix */
   nr_rotzf(TRZ, rot_z);             /* create the rotate Z matrix */
   
   nr_multf(TRY,1,4,1,4,  TRX,1,4,1,4,  T1); /* apply rx and ry */
   nr_multf(TRZ,1,4,1,4,  T1,1,4,1,4,   xmat); /* apply rz */


   free_matrix(TRX,1,4,1,4);
   free_matrix(TRY,1,4,1,4);
   free_matrix(TRZ,1,4,1,4);
   free_matrix(T1 ,1,4,1,4);
   free_matrix(T2 ,1,4,1,4);

}




/* ----------------------------- MNI Header -----------------------------------
@NAME       : build_transformation_matrix
@INPUT      : center, translations, scales, rotations
@OUTPUT     : *lt->mat - a linear transformation matrix
@RETURNS    : nothing
@DESCRIPTION: mat = (T)(C)(SH)(S)(R)(-C)
               the matrix is to be  PREmultiplied with a column vector (mat*colvec)
	       when used in the application
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Thu Jun  3 09:37:56 EST 1993 lc
@MODIFIED   : 
---------------------------------------------------------------------------- */

public void build_transformation_matrix(Transform *trans,
					double *center,
					double *translations,
					double *scales,
					double *shears,
					double *rotations)
{
  
  float
    **T,
    **SH,
    **S,
    **R,
    **C,
    **T1,
    **T2,
    **T3,
    **T4;
  int
    i,j;
  
  T  = matrix(1,4,1,4);
  SH = matrix(1,4,1,4);
  S  = matrix(1,4,1,4);
  R  = matrix(1,4,1,4);
  C  = matrix(1,4,1,4);
  T1 = matrix(1,4,1,4);
  T2 = matrix(1,4,1,4);
  T3 = matrix(1,4,1,4);
  T4 = matrix(1,4,1,4);
  
				             /* mat = (T)(C)(SH)(S)(R)(-C) */

  nr_identf(T,1,4,1,4);	             /* make (T)(C) */
  for_less( i, 0, 3 ) {
    T[1+i][4] = translations[i] - center[i];		
  }
                                /* make rotation matix */
  make_rots(R,
	    (float)(rotations[0]),
	    (float)(rotations[1]),
	    (float)(rotations[2])); 

                                /* make shear rotation matrix */
  make_rots(SH,
	    (float)(shears[0]),
	    (float)(shears[1]),      
	    (float)(shears[2])); 

                                /* make scaling matrix */
  nr_identf(S,1,4,1,4);	           
  for_less( i, 0, 3 ) {
    S[1+i][1+i] = scales[i];
  }

  nr_identf(C,1,4,1,4);      /* make center          */
  for_less( i, 0, 3 ) {
    C[1+i][4] = center[i];		
  }

  nr_multf(C,1,4,1,4,    SH  ,1,4,1,4, T1 );  
  nr_multf(T1,1,4,1,4,   S ,1,4,1,4, T2 );  
  nr_multf(T2,1,4,1,4,   R   ,1,4,1,4, T3 );  
  nr_multf(T3,1,4,1,4,   T,1,4,1,4, T4 );  

  for_less(i,0,4)
    for_less(j,0,4)
      Transform_elem(*trans, i, j ) = T4[i+1][j+1];

  free_matrix(T    ,1,4,1,4);
  free_matrix(SH   ,1,4,1,4);
  free_matrix(S    ,1,4,1,4);
  free_matrix(R    ,1,4,1,4);
  free_matrix(C ,1,4,1,4);
  free_matrix(T1   ,1,4,1,4);
  free_matrix(T2   ,1,4,1,4);
  free_matrix(T3   ,1,4,1,4);
  free_matrix(T4   ,1,4,1,4);
}



public void  build_homogeneous_from_parameters(float **calc_transformation,
						int ndim, 
						float *centre, 
						float *translation,
						float *scales,
						float *shears,
						float *angles)
{
  double cent[3];
  double trans[3];
  double sc[3];
  double sh[3];
  double ang[3];
  Transform mat;

  int i,j;
  
  for(i=0; i<ndim; i++) {
    cent[i]  = (double)centre[i+1];
    trans[i] = (double)translation[i+1];
    sc[i]    = (double)scales[i+1];
    sh[i]    = (double)shears[i+1];
    ang[i]   = (double)angles[i+1];
  }

  build_transformation_matrix(&mat,cent, trans, sc, sh, ang);

  for (i=0; i<ndim; i++) 
    for (j=0; j<=ndim; j++)
      calc_transformation[j+1][i+1] = Transform_elem(mat,i,j);

  for (i=0; i<ndim; i++) 
    calc_transformation[i+1][ndim+1] = 0;

  calc_transformation[ndim+1][ndim+1] = 1.0;

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : build_inverse_transformation_matrix
@INPUT      : center, translations, scales, rotations
@OUTPUT     : the inverse linear transformation matrix of mat:
                since mat = (T)(C)(SH)(S)(R)(-C), then

		invmat = (C)(inv(r))(inv(S))(inv(SH))(-C)(-T)

@RETURNS    : nothing
@DESCRIPTION: 
               the matrix is to be  PREmultiplied with a vector (mat*vec)
	       when used in the application
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Tue Jun 15 16:45:35 EST 1993 LC
@MODIFIED   : 
---------------------------------------------------------------------------- */
public void build_inverse_transformation_matrix(Transform *trans,
						double *center,
						double *translations,
						double *scales,
						double *shears,
						double *rotations)
{
  float
    **T,
    **SH,
    **S,
    **R,
    **C,
    **T1,
    **T2,
    **T3,
    **T4;
  int
    i,j;
  
  T      = matrix(1,4,1,4);
  SH     = matrix(1,4,1,4);
  S      = matrix(1,4,1,4);
  R      = matrix(1,4,1,4);
  C   = matrix(1,4,1,4);
  T1     = matrix(1,4,1,4);
  T2     = matrix(1,4,1,4);
  T3     = matrix(1,4,1,4);
  T4     = matrix(1,4,1,4);
  
				/* invmat = (C)(inv(r))(inv(S))(inv(SH))(-C)(-T)
				   mat = (T)(C)(SH)(S)(R)(-C) */

  nr_identf(T,1,4,1,4);	             /* make (-T)(-C) */
  for_less( i, 0, 3 ) {
    T[1+i][4] = -translations[i] + center[i];		
  }
                                /* make rotation matix */
  make_rots(T1,
	    (float)(rotations[0]),
	    (float)(rotations[1]),
	    (float)(rotations[2])); 
  transpose(4,4,T1,R);

                                /* make shear rotation matrix */
  make_rots(T1,
	    (float)(shears[0]),
	    (float)(shears[1]),      
	    (float)(shears[2])); 
  transpose(4,4,T1,SH);

                                /* make scaling matrix */
  nr_identf(S,1,4,1,4);	           
  for_less( i, 0, 3 ) {
    if (scales[i] != 0.0)
      S[1+i][1+i] = 1/scales[i];
    else
      S[1+i][1+i] = 1.0;
  }

  nr_identf(C,1,4,1,4);      /* make center          */
  for_less( i, 0, 3 ) {
    C[1+i][4] = -center[i];		
  }

  nr_multf(T,1,4,1,4,  R ,1,4,1,4, T1 );  
  nr_multf(T1,1,4,1,4, S ,1,4,1,4, T2 );  
  nr_multf(T2,1,4,1,4, SH,1,4,1,4, T3 );  
  nr_multf(T3,1,4,1,4, C ,1,4,1,4, T4 );  

  for_less(i,0,4)
    for_less(j,0,4)
      Transform_elem(*trans, i, j ) = T4[i+1][j+1];

  free_matrix(T    ,1,4,1,4);
  free_matrix(SH   ,1,4,1,4);
  free_matrix(S    ,1,4,1,4);
  free_matrix(R    ,1,4,1,4);
  free_matrix(C ,1,4,1,4);
  free_matrix(T1   ,1,4,1,4);
  free_matrix(T2   ,1,4,1,4);
  free_matrix(T3   ,1,4,1,4);
  free_matrix(T4   ,1,4,1,4);
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : extract_parameters_from_matrix
@INPUT      : trans    - a linear transformation matrix structure
              center   - an array of the desired center of rotation and scaling.
@OUTPUT     : translations, scales, rotations
@RETURNS    : nothing
@DESCRIPTION: mat = (C)(SH)(S)(R)(-C)(T)
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Thu Jun  3 09:37:56 EST 1993 lc
@MODIFIED   : 
---------------------------------------------------------------------------- */

public BOOLEAN extract_parameters_from_matrix(Transform *trans,
					      double *center,
					      double *translations,
					      double *scales,
					      double *shears,
					      double *rotations)
{
  int 
    i,j;

  float 
    magnitude,
    **center_of_rotation,
    **result,
    **unit_vec,
    *ang,*tmp,
    **xmat,**T,**Tinv,**C,**Sinv,**R,**SR,**SRinv,**Cinv,**TMP1,**TMP2;

  xmat  = matrix(1,4,1,4); nr_identf(xmat ,1,4,1,4);
  TMP1  = matrix(1,4,1,4); nr_identf(TMP1 ,1,4,1,4);
  TMP2  = matrix(1,4,1,4); nr_identf(TMP2 ,1,4,1,4);
  Cinv  = matrix(1,4,1,4); nr_identf(Cinv ,1,4,1,4);
  SR    = matrix(1,4,1,4); nr_identf(SR   ,1,4,1,4);
  SRinv = matrix(1,4,1,4); nr_identf(SRinv,1,4,1,4);
  Sinv  = matrix(1,4,1,4); nr_identf(Sinv ,1,4,1,4); 
  T     = matrix(1,4,1,4); nr_identf(T    ,1,4,1,4);
  Tinv  = matrix(1,4,1,4); nr_identf(Tinv ,1,4,1,4);
  C     = matrix(1,4,1,4); nr_identf(C    ,1,4,1,4);
  R     = matrix(1,4,1,4); nr_identf(R    ,1,4,1,4);

  center_of_rotation = matrix(1,4,1,1);	/* make column vectors */
  result             = matrix(1,4,1,1);
  unit_vec           = matrix(1,4,1,1);

  tmp = vector(1,3);
  ang = vector(1,3);


  for_inclusive( i, 0, 3)	/* copy the input matrix */
    for_inclusive( j, 0, 3)
      xmat[i+1][j+1] = (float)Transform_elem(*trans,i,j);

  

  /* -------------DETERMINE THE TRANSLATION FIRST! ---------  */

				/* see where the center of rotation is displaced... */

  FILL_NR_COLVEC( center_of_rotation, center[0], center[1], center[2] );


  invertmatrix(4, xmat, TMP1);	/* get inverse of the matrix */

  raw_matrix_multiply( 4, 4, 1, TMP1, center_of_rotation, result);

  SUB_NR_COLVEC( result, result, center_of_rotation );

  for_inclusive( i, 0, 2) 
    translations[i] = -result[i+1][1];

  /* -------------NOW GET THE SCALING VALUES! ----------------- */

  for_inclusive( i, 0, 2) 
    tmp[i+1] = -translations[i];
  translation_to_homogeneous(3, tmp, Tinv); 

  for_inclusive( i, 0, 2) 
    tmp[i+1] = center[i];
  translation_to_homogeneous(3, tmp, C); 
  for_inclusive( i, 0, 2) 
    tmp[i+1] = -center[i];
  translation_to_homogeneous(3, tmp, Cinv); 


  matrix_multiply(4,4,4, Cinv, xmat, TMP1);    /* get scaling*rotation matrix */


  matrix_multiply(4,4,4, TMP1, C, TMP1);


  matrix_multiply(4,4,4, TMP1, Tinv,    SR);

  invertmatrix(4, SR, SRinv);	/* get inverse of scaling*rotation */

				/* find each scale by mapping a unit vector backwards,
				   and finding the magnitude of the result. */
  FILL_NR_COLVEC( unit_vec, 1.0, 0.0, 0.0 );
  raw_matrix_multiply( 4, 4, 1, SRinv, unit_vec, result);
  magnitude = MAG_NR_COLVEC( result );
  if (magnitude != 0.0) {
    scales[0] = 1/magnitude;
    Sinv[1][1] = magnitude;
  }
  else {
    scales[0] = 1.0;
    Sinv[1][1] = 1.0;
  }

  FILL_NR_COLVEC( unit_vec, 0.0, 1.0, 0.0 );
  raw_matrix_multiply( 4, 4, 1, SRinv, unit_vec, result);
  magnitude = MAG_NR_COLVEC( result );
  if (magnitude != 0.0) {
    scales[1] = 1/magnitude;
    Sinv[2][2] = magnitude;
  }
  else {
    scales[1]  = 1.0;
    Sinv[2][2] = 1.0;
  }

  FILL_NR_COLVEC( unit_vec, 0.0, 0.0, 1.0 );
  raw_matrix_multiply( 4, 4, 1, SRinv, unit_vec, result);
  magnitude = MAG_NR_COLVEC( result );
  if (magnitude != 0.0) {
    scales[2] = 1/magnitude;
    Sinv[3][3] = magnitude;
  }
  else {
    scales[2] = 1.0;
    Sinv[3][3] = 1.0;
  }

  /* ------------NOW GET THE ROTATION ANGLES!----- */

				/* extract rotation matrix */
  matrix_multiply(4,4, 4, Sinv, SR,   R);

				/* get rotation angles */
  if (!rotmat_to_ang(R, ang)) {
    (void)fprintf(stderr,"Cannot convert R to radians!");
    printmatrix(3,3,R);
    return(FALSE);
  }

  for_inclusive( i, 0, 2 )
    rotations[i] = ang[i+1];

  free_matrix(xmat,1,4,1,4);
  free_matrix(TMP1,1,4,1,4);
  free_matrix(TMP2,1,4,1,4);
  free_matrix(Cinv,1,4,1,4);
  free_matrix(SR  ,1,4,1,4);
  free_matrix(SRinv,1,4,1,4);
  free_matrix(Sinv,1,4,1,4);
  free_matrix(T   ,1,4,1,4);
  free_matrix(Tinv,1,4,1,4);
  free_matrix(C   ,1,4,1,4);
  free_matrix(R   ,1,4,1,4);
  
  free_matrix(center_of_rotation,1,4,1,1);	/* make column vectors */
  free_matrix(result            ,1,4,1,1);
  free_matrix(unit_vec          ,1,4,1,1);

  free_vector(ang, 1, 3);
  free_vector(tmp, 1, 3);

  return(TRUE);
}

















