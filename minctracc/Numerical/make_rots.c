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
@MODIFIED   : Revision 96.7  2006-11-29 09:09:33  rotor
@MODIFIED   :  * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   : Revision 96.6  2005/07/20 20:45:49  rotor
@MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :     * Removed old VOLUME_IO cruft #defines
@MODIFIED   :     * Fixed up all Makefile.am's in subdirs
@MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   : Revision 96.5  2004/02/12 05:54:27  rotor
@MODIFIED   :  * removed /static defs
@MODIFIED   :
@MODIFIED   : Revision 96.4  2002/11/20 21:38:49  lenezet
@MODIFIED   :
@MODIFIED   : Fix the code to take in consideration the direction cosines especially in the grid transform.
@MODIFIED   : Add an option to choose the maximum expected deformation magnitude.
@MODIFIED   :
@MODIFIED   : Revision 96.3  2002/08/14 19:54:26  lenezet
@MODIFIED   :  quaternion option added for the rotation
@MODIFIED   :
@MODIFIED   : Revision 96.2  2002/03/26 14:15:40  stever
@MODIFIED   : Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   : Revision 96.1  1997/11/03 19:59:49  louis
@MODIFIED   : - now include volume_io/internal_volume_io.h instead of volume_io.h
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
 * Revision 1.12  1996/08/12  14:15:45  louis
 * Pre-release
 *
 * Revision 1.11  1995/09/11  12:37:16  collins
 * All refs to numerical recipes routines have been replaced.
 * this is an updated working version - corresponds to mni_reg-0.1g
 * \
 *
 * Revision 1.10  1995/02/22  08:56:06  collins
 * Montreal Neurological Institute version.
 * compiled and working on SGI.  this is before any changes for SPARC/
 * Solaris.
 *
 * Revision 1.9  94/04/26  12:54:26  louis
 * updated with new versions of make_rots, extract2_parameters_from_matrix 
 * that include proper interpretation of skew.
 * 
 * Revision 1.8  94/04/06  11:48:39  louis
 * working linted version of linear + non-linear registration based on Lvv
 * operator working in 3D
 * 
 * Revision 1.7  94/02/21  16:35:42  louis
 * version before feb 22 changes
 * 
 * Revision 1.6  93/11/15  16:27:02  louis
 * working version, with new library, with RCS revision stuff,
 * before deformations included
 * 
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /static-cvsroot/registration/mni_autoreg/minctracc/Numerical/make_rots.c,v 96.7 2006-11-29 09:09:33 rotor Exp $";
#endif

#include "config.h"

#include <volume_io.h>
#include "matrix_basics.h"
#include "rotmat_to_ang.h"
#include "local_macros.h"

#define  FILL_NR_COLVEC( vector, x, y, z ) \
            { \
                vector[1][1] = (x); \
                vector[2][1] = (y); \
                vector[3][1] = (z); \
                vector[4][1] = 1.0; \
            }

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

#define  DOT_NR_COLVEC( vector, v1, v2 ) \
            { \
                vector[1][1] = v1[1][1]*v2[1][1]; \
                vector[2][1] = v1[2][1]*v2[2][1]; \
                vector[3][1] = v1[3][1]*v2[3][1]; \
                vector[4][1] = 1.0; \
            }

#define  SCALAR_MULT_NR_COLVEC( vector, v1, sc ) \
            { \
                vector[1][1] = v1[1][1]*sc; \
                vector[2][1] = v1[2][1]*sc; \
                vector[3][1] = v1[3][1]*sc; \
                vector[4][1] = 1.0; \
            }

#define  DOTSUM_NR_COLVEC( v1, v2 ) \
                (v1[1][1] * v2[1][1] + \
                v1[2][1] * v2[2][1] + \
                v1[3][1] * v2[3][1]) 

#define  MAG_NR_COLVEC( v1 ) \
            ( sqrt( v1[1][1] * v1[1][1] + \
                v1[2][1] * v1[2][1] + \
                v1[3][1] * v1[3][1] ) )


void build_rotmatrix(float **m, double *quat);
void extract_quaternions(float **m, double *quat);

/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_rots
@INPUT      : rot_x, rot_y, rot_z - three rotation angles, in radians.
@OUTPUT     : xmat, a zero offset matrix for homogeous transformations
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
      **T1;
   
   VIO_ALLOC2D(TRX  ,5,5);
   VIO_ALLOC2D(TRY  ,5,5);
   VIO_ALLOC2D(TRZ  ,5,5);
   VIO_ALLOC2D(T1   ,5,5);

   nr_rotxf(TRX, rot_x);             /* create the rotate X matrix */
   nr_rotyf(TRY, rot_y);             /* create the rotate Y matrix */
   nr_rotzf(TRZ, rot_z);             /* create the rotate Z matrix */
   
   nr_multf(TRY,1,4,1,4,  TRX,1,4,1,4,  T1); /* apply rx and ry */
   nr_multf(TRZ,1,4,1,4,  T1,1,4,1,4,   xmat); /* apply rz */


   VIO_FREE2D(TRX);
   VIO_FREE2D(TRY);
   VIO_FREE2D(TRZ);
   VIO_FREE2D(T1 );

}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_shears
@INPUT      : shear - an array of six shear parameters.
@OUTPUT     : xmat, a zero offset matrix for homogeous transformations
@RETURNS    : nothing
@DESCRIPTION: to be applied by premultiplication, ie rot*vec = newvec
@METHOD     : 
                xmat = [ 1 a b 0
                         c 1 d 0
                         e f 1 0
                         0 0 0 1 ];

                where shear[] = [a b c d e f]     

@CREATED    : Sat Apr 16 10:44:26 EST 1994
@MODIFIED   : 
---------------------------------------------------------------------------- */

void   make_shears(float **xmat,                                         
                 double *shears)
{

    nr_identf(xmat,1,4,1,4);        /* start with identity */
    xmat[2][1] = shears[0];
    xmat[3][1] = shears[1];
    xmat[3][2] = shears[2];
  
}




/* ----------------------------- MNI Header -----------------------------------
@NAME       : build_transformation_matrix
@INPUT      : center, translations, scales, rotations
@OUTPUT     : *lt->mat - a linear transformation matrix
@RETURNS    : nothing
@DESCRIPTION: mat = (T)(C)(S)(SH)(R)(-C)
               the matrix is to be  PREmultiplied with a column vector (mat*colvec)
               when used in the application
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Thu Jun  3 09:37:56 EST 1993 lc
@MODIFIED   : 
---------------------------------------------------------------------------- */

void build_transformation_matrix(VIO_Transform *trans,
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
  
  VIO_ALLOC2D(T  ,5,5);
  VIO_ALLOC2D(SH ,5,5);
  VIO_ALLOC2D(S  ,5,5);
  VIO_ALLOC2D(R  ,5,5);
  VIO_ALLOC2D(C  ,5,5);
  VIO_ALLOC2D(T1 ,5,5);
  VIO_ALLOC2D(T2 ,5,5);
  VIO_ALLOC2D(T3 ,5,5);
  VIO_ALLOC2D(T4 ,5,5);
  
                                             /* mat = (T)(C)(SH)(S)(R)(-C) */

  nr_identf(T,1,4,1,4);                     /* make (T)(C) */
  for(i=0; i<3; i++) {
    T[1+i][4] = translations[i] + center[i];                
  }
                                /* make rotation matix */
  make_rots(R,
            (float)(rotations[0]),
            (float)(rotations[1]),
            (float)(rotations[2])); 

                                /* make shear rotation matrix */
  make_shears(SH, shears);

                                /* make scaling matrix */
  nr_identf(S,1,4,1,4);                   
  for(i=0; i<3; i++) {
    S[1+i][1+i] = scales[i];
  }

  nr_identf(C,1,4,1,4);      /* make center          */
  for(i=0; i<3; i++) {
    C[1+i][4] = -center[i];                
  }

  nr_multf(T, 1,4,1,4, S  ,1,4,1,4, T1 );  
  nr_multf(T1,1,4,1,4, SH ,1,4,1,4, T2 );  
  nr_multf(T2,1,4,1,4, R  ,1,4,1,4, T3 );  
  nr_multf(T3,1,4,1,4, C  ,1,4,1,4, T4 );  

  for(i=0; i<4; i++)
    for(j=0; j<4; j++)
      Transform_elem(*trans, i, j ) = T4[i+1][j+1];

  VIO_FREE2D(T    );
  VIO_FREE2D(SH   );
  VIO_FREE2D(S    );
  VIO_FREE2D(R    );
  VIO_FREE2D(C );
  VIO_FREE2D(T1   );
  VIO_FREE2D(T2   );
  VIO_FREE2D(T3   );
  VIO_FREE2D(T4   );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : build_transformation_matrix_quater
@INPUT      : center, translations, scales, quaternions
@OUTPUT     : *lt->mat - a linear transformation matrix
@RETURNS    : nothing
@DESCRIPTION: mat = (T)(C)(S)(SH)(R)(-C)
               the matrix is to be  PREmultiplied with a column vector (mat*colvec)
               when used in the application
               same as build_transformation_matrix but with quaternions
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Thr Apr 18 10:45:56 EST 2002 pln
@MODIFIED   : 
---------------------------------------------------------------------------- */

void build_transformation_matrix_quater(VIO_Transform *trans,
                                               double *center,
                                               double *translations,
                                               double *scales,
                                               double *shears,
                                               double *quaternions)
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

  double normal;


  int
    i,j;
  
  VIO_ALLOC2D(T  ,5,5);
  VIO_ALLOC2D(SH ,5,5);
  VIO_ALLOC2D(S  ,5,5);
  VIO_ALLOC2D(R  ,5,5);
  VIO_ALLOC2D(C  ,5,5);
  VIO_ALLOC2D(T1 ,5,5);
  VIO_ALLOC2D(T2 ,5,5);
  VIO_ALLOC2D(T3 ,5,5);
  VIO_ALLOC2D(T4 ,5,5);
  
 
  


  
  normal=(quaternions[0]*quaternions[0] + quaternions[1]*quaternions[1] + quaternions[2]*quaternions[2] + quaternions[3]*quaternions[3]);
  if (normal>1){
   for(i = 0; i < 4; i++){
      quaternions[i] /= normal;
   }} 
   
   

                        /* mat = (T)(C)(SH)(S)(R)(-C) */

  nr_identf(T,1,4,1,4);



                                             /* make (T)(C) */
  for(i=0; i<3; i++) {
    T[1+i][4] = translations[i] + center[i];                
  }
   
  

  build_rotmatrix(R,quaternions); /* make rotation matrix from quaternions */
  

                                /* make shear rotation matrix */
  make_shears(SH, shears);

                                /* make scaling matrix */
  nr_identf(S,1,4,1,4);                   
  for(i=0; i<3; i++) {
    S[1+i][1+i] = scales[i];
  }


  nr_identf(C,1,4,1,4);      /* make center          */
  for(i=0; i<3; i++) {
    C[1+i][4] = -center[i];                
  }

  nr_multf(T, 1,4,1,4, S  ,1,4,1,4, T1 );  
  nr_multf(T1,1,4,1,4, SH ,1,4,1,4, T2 );  
  nr_multf(T2,1,4,1,4, R  ,1,4,1,4, T3 );  
  nr_multf(T3,1,4,1,4, C  ,1,4,1,4, T4 );  

  for(i=0; i<4; i++)
    for(j=0; j<4; j++)
      Transform_elem(*trans, i, j ) = T4[i+1][j+1];

  VIO_FREE2D(T    );
  VIO_FREE2D(SH   );
  VIO_FREE2D(S    );
  VIO_FREE2D(R    );
  VIO_FREE2D(C );
  VIO_FREE2D(T1   );
  VIO_FREE2D(T2   );
  VIO_FREE2D(T3   );
  VIO_FREE2D(T4   );
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
void build_inverse_transformation_matrix(VIO_Transform *trans,
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
  
  VIO_ALLOC2D(T   ,5,5);
  VIO_ALLOC2D(SH  ,5,5);
  VIO_ALLOC2D(S   ,5,5);
  VIO_ALLOC2D(R   ,5,5);
  VIO_ALLOC2D(C   ,5,5);
  VIO_ALLOC2D(T1  ,5,5);
  VIO_ALLOC2D(T2  ,5,5);
  VIO_ALLOC2D(T3  ,5,5);
  VIO_ALLOC2D(T4  ,5,5);
  
                                /* invmat = (C)(inv(r))(inv(S))(inv(SH))(-C)(-T)
                                   mat = (T)(C)(SH)(S)(R)(-C) */

  nr_identf(T,1,4,1,4);                     /* make (-T)(-C) */
  for(i=0; i<3; i++) {
    T[1+i][4] = -translations[i] - center[i];                
  }

                                /* make rotation matix */
  make_rots(T1,
            (float)(rotations[0]),
            (float)(rotations[1]),
            (float)(rotations[2])); 


  transpose(4,4,T1,R);
  

  make_shears(T1,shears);        /* make shear rotation matrix */
  invertmatrix(4, T1, SH);        /* get inverse of the matrix */

                                /* make scaling matrix */
  nr_identf(S,1,4,1,4);                   
  for(i=0; i<3; i++) {
    if (scales[i] != 0.0)
      S[1+i][1+i] = 1/scales[i];
    else
      S[1+i][1+i] = 1.0;
  }

  nr_identf(C,1,4,1,4);      /* make center          */
  for(i=0; i<3; i++) {
    C[1+i][4] = center[i];                
  }

  nr_multf(C,1,4,1,4,  R ,1,4,1,4, T1 );  
  nr_multf(T1,1,4,1,4, SH,1,4,1,4, T2 );  
  nr_multf(T2,1,4,1,4, S ,1,4,1,4, T3 );  
  nr_multf(T3,1,4,1,4, T ,1,4,1,4, T4 );  

  for(i=0; i<4; i++)
    for(j=0; j<4; j++)
      Transform_elem(*trans, i, j ) = T4[i+1][j+1];

  VIO_FREE2D(T    );
  VIO_FREE2D(SH   );
  VIO_FREE2D(S    );
  VIO_FREE2D(R    );
  VIO_FREE2D(C );
  VIO_FREE2D(T1   );
  VIO_FREE2D(T2   );
  VIO_FREE2D(T3   );
  VIO_FREE2D(T4   );
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : build_inverse_transformation_matrix_quater
@INPUT      : center, translations, scales, quaternions
@OUTPUT     : the inverse linear transformation matrix of mat:
                since mat = (T)(C)(SH)(S)(R)(-C), then

                invmat = (C)(inv(r))(inv(S))(inv(SH))(-C)(-T)

@RETURNS    : nothing
@DESCRIPTION: 
               the matrix is to be  PREmultiplied with a vector (mat*vec)
               when used in the application
               same as build_inverse_transformation_matrix but with quaternions
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Thr Apr 18 10:45:56 EST 2002 pln
@MODIFIED   : 
---------------------------------------------------------------------------- */
void build_inverse_transformation_matrix_quater(VIO_Transform *trans,
                                                       double *center,
                                                       double *translations,
                                                       double *scales,
                                                       double *shears,
                                                       double *quaternions)
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
  
  VIO_ALLOC2D(T   ,5,5);
  VIO_ALLOC2D(SH  ,5,5);
  VIO_ALLOC2D(S   ,5,5);
  VIO_ALLOC2D(R   ,5,5);
  VIO_ALLOC2D(C   ,5,5);
  VIO_ALLOC2D(T1  ,5,5);
  VIO_ALLOC2D(T2  ,5,5);
  VIO_ALLOC2D(T3  ,5,5);
  VIO_ALLOC2D(T4  ,5,5);
  
                                /* invmat = (C)(inv(r))(inv(S))(inv(SH))(-C)(-T)
                                   mat = (T)(C)(SH)(S)(R)(-C) */

  nr_identf(T,1,4,1,4);                     /* make (-T)(-C) */
  for(i=0; i<3; i++) {
    T[1+i][4] = -translations[i] - center[i];                
  }

  

  build_rotmatrix(T1,quaternions); /* make rotation matrix from quaternions */
  transpose(4,4,T1,R);
  

  make_shears(T1,shears);        /* make shear rotation matrix */
  invertmatrix(4, T1, SH);        /* get inverse of the matrix */

                                /* make scaling matrix */
  nr_identf(S,1,4,1,4);                   
  for(i=0; i<3; i++) {
    if (scales[i] != 0.0)
      S[1+i][1+i] = 1/scales[i];
    else
      S[1+i][1+i] = 1.0;
  }

  nr_identf(C,1,4,1,4);      /* make center          */
  for(i=0; i<3; i++) {
    C[1+i][4] = center[i];                
  }

  nr_multf(C,1,4,1,4,  R ,1,4,1,4, T1 );  
  nr_multf(T1,1,4,1,4, SH,1,4,1,4, T2 );  
  nr_multf(T2,1,4,1,4, S ,1,4,1,4, T3 );  
  nr_multf(T3,1,4,1,4, T ,1,4,1,4, T4 );  

  for(i=0; i<4; i++)
    for(j=0; j<4; j++)
      Transform_elem(*trans, i, j ) = T4[i+1][j+1];

  VIO_FREE2D(T    );
  VIO_FREE2D(SH   );
  VIO_FREE2D(S    );
  VIO_FREE2D(R    );
  VIO_FREE2D(C );
  VIO_FREE2D(T1   );
  VIO_FREE2D(T2   );
  VIO_FREE2D(T3   );
  VIO_FREE2D(T4   );
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
@MODIFIED   : Sun Apr 17 09:54:14 EST 1994 - tried to extract shear parameters

      if det(ROT) != 1, then the ROT matrix is not a pure rotation matrix.
      I will find the shear matrix required, by building a rotation matrix R1
      from the extracted rotation parameters (rx,ry and rz), and multiply ROT by 
      inv(R1).
---------------------------------------------------------------------------- */

VIO_BOOL extract_parameters_from_matrix(VIO_Transform *trans,
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

  VIO_ALLOC2D(xmat  ,5,5); nr_identf(xmat ,1,4,1,4);
  VIO_ALLOC2D(TMP1  ,5,5); nr_identf(TMP1 ,1,4,1,4);
  VIO_ALLOC2D(TMP2  ,5,5); nr_identf(TMP2 ,1,4,1,4);
  VIO_ALLOC2D(Cinv  ,5,5); nr_identf(Cinv ,1,4,1,4);
  VIO_ALLOC2D(SR    ,5,5); nr_identf(SR   ,1,4,1,4);
  VIO_ALLOC2D(SRinv ,5,5); nr_identf(SRinv,1,4,1,4);
  VIO_ALLOC2D(Sinv  ,5,5); nr_identf(Sinv ,1,4,1,4); 
  VIO_ALLOC2D(T     ,5,5); nr_identf(T    ,1,4,1,4);
  VIO_ALLOC2D(Tinv  ,5,5); nr_identf(Tinv ,1,4,1,4);
  VIO_ALLOC2D(C     ,5,5); nr_identf(C    ,1,4,1,4);
  VIO_ALLOC2D(R     ,5,5); nr_identf(R    ,1,4,1,4);

  VIO_ALLOC2D(center_of_rotation ,5,5);        /* make column vectors */
  VIO_ALLOC2D(result             ,5,5);
  VIO_ALLOC2D(unit_vec           ,5,5);

  ALLOC(tmp ,4);
  ALLOC(ang ,4);


  for(i=0; i<=3; i++)        /* copy the input matrix */
    for(j=0; j<=3; j++)
      xmat[i+1][j+1] = (float)Transform_elem(*trans,i,j);

  

  /* -------------DETERMINE THE TRANSLATION FIRST! ---------  */

                                /* see where the center of rotation is displaced... */

  FILL_NR_COLVEC( center_of_rotation, center[0], center[1], center[2] );


  invertmatrix(4, xmat, TMP1);        /* get inverse of the matrix */

  matrix_multiply( 4, 4, 1, xmat, center_of_rotation, result); /* was TMP! in place of xmat */

  SUB_NR_COLVEC( result, result, center_of_rotation );

  for(i=0; i<=2; i++) 
    translations[i] = result[i+1][1];

  /* -------------NOW GET THE SCALING VALUES! ----------------- */

  for(i=0; i<=2; i++) 
    tmp[i+1] = -translations[i];
  translation_to_homogeneous(3, tmp, Tinv); 

  for(i=0; i<=2; i++) 
    tmp[i+1] = center[i];
  translation_to_homogeneous(3, tmp, C); 
  for(i=0; i<=2; i++) 
    tmp[i+1] = -center[i];
  translation_to_homogeneous(3, tmp, Cinv); 


  matrix_multiply(4,4,4, xmat, C, TMP1);    /* get scaling*rotation matrix */


  matrix_multiply(4,4,4, Tinv, TMP1, TMP1);


  matrix_multiply(4,4,4, Cinv, TMP1,    SR);

  invertmatrix(4, SR, SRinv);        /* get inverse of scaling*rotation */

                                /* find each scale by mapping a unit vector backwards,
                                   and finding the magnitude of the result. */
  FILL_NR_COLVEC( unit_vec, 1.0, 0.0, 0.0 );
  matrix_multiply( 4, 4, 1, SRinv, unit_vec, result);
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
  matrix_multiply( 4, 4, 1, SRinv, unit_vec, result);
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
  matrix_multiply( 4, 4, 1, SRinv, unit_vec, result);
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

  for(i=0; i<=2; i++)
    rotations[i] = ang[i+1];


  VIO_FREE2D(xmat);
  VIO_FREE2D(TMP1);
  VIO_FREE2D(TMP2);
  VIO_FREE2D(Cinv);
  VIO_FREE2D(SR  );
  VIO_FREE2D(SRinv);
  VIO_FREE2D(Sinv);
  VIO_FREE2D(T   );
  VIO_FREE2D(Tinv);
  VIO_FREE2D(C   );
  VIO_FREE2D(R   );
  
  VIO_FREE2D(center_of_rotation);
  VIO_FREE2D(result            );
  VIO_FREE2D(unit_vec          );

  FREE(ang);
  FREE(tmp);

  return(TRUE);
}




/*
   function getparams will get the trnsform  parameters from a 
   transformation matrix 'tran' (that has already had the translation 
   componants removed).
   in description below, I assume that tran is a forward xform, from 
   native space to talairach space.
   I assume that trans is a 4x4 homogeneous matrix, with the 
   principal axis stored in the upper left 3x3 matrix.  The first col
   of tran represents is the coordinate of (1,0,0)' mapped through the
   transformation.  (this means  vec2 = tran * vec1).
  
   trans = [scale][shear][rot]
         = [scale][shear][rz][ry][rx];

   the shear matrix is constrained to be of the form:
     shear = [1 0 0 0
              f 1 0 0
              g h 1 0
              0 0 0 1];
     where f,g,h can take on any value.

   the scale matrix is constrained to be of the form:
     scale = [sx 0  0  0
              0  sy 0  0
              0  0  sz 0
              0  0  0  1];
     where sx,sy,sz must be positive.
   
  all rotations are assumed to be in the range -pi/2..pi/2

  the rotation angles are returned as radians and are applied 
  counter clockwise, when looking down the axis (from the positive end
  towards the origin).

  trans is assumed to be invertible.

  i assume a coordinate system:
             ^ y
             |
             |
             |
             |_______> x
            /
           /
          /z  (towards the viewer).

  
  
  procedure: 
          start with t = inv(tran) (after removing translation from tran )
  
          t = inv(r)*inv(sh)*inv(sc)
                  t maps the talairach space unit vectors into native
                space, which means that the columns of T are the
                direction cosines of these vectors x,y,z

   (1)  the length of the vectors x,y,z give the inverse scaling
        parameters:
             sx = 1/norm(x); sy = 1/norm(y); sz = 1/norm(z);

   (2)  with the constraints on the form of sh above, inv(sh) has the
        same form.  let inv(sh) be as above in terms of a,b and c.
          inv(sh) = [1 0 0 0; a 1 0 0; b c 1 0; 0 0 0 1];

        for c: project y onto z and normalize:

                  /     |y.z|^2     \(1/2)
             c = <  ---------------  >
                  \ |y|^2 - |y.z|^2 /

        for b: project x onto z and normalize

                  /        |x.z|^2        \(1/2)
             b = <  ---------------------  >
                  \ |x|^2 - |x.z|^2 - a^2 /
     
          where a is the projection of x onto the coordinate sys Y axis.

        for a: project x onto z and normalize
           a is taken from (b) above, and normalized... see below

  (3) rots are returned by getrots by giving the input transformation:
        rot_mat = [inv(sh)][inv(sc)][trans]

  (4) once completed, the parameters of sx,sy,sz and a,b,c are
      adjusted so that they maintain the matrix contraints above.


*/
VIO_BOOL extract2_parameters_from_matrix(VIO_Transform *trans,
                                               double *center,
                                               double *translations,
                                               double *scales,
                                               double *shears,
                                               double *rotations)
{
  int 
    i,j;

  float 
    n1,n2,
    magnitude, magz, magx, magy, ai,bi,ci,scalar,a1,
    **center_of_rotation,
    **result,
    **unit_vec,
    *ang,*tmp,**x,**y,**z, **nz, **y_on_z, **ortho_y,
    **xmat,**T,**Tinv,**C,**Sinv,
    **R,**SR,**SRinv,**Cinv,**TMP1,**TMP2;

  VIO_ALLOC2D(xmat  ,5,5); nr_identf(xmat ,1,4,1,4);
  VIO_ALLOC2D(TMP1  ,5,5); nr_identf(TMP1 ,1,4,1,4);
  VIO_ALLOC2D(TMP2  ,5,5); nr_identf(TMP2 ,1,4,1,4);
  VIO_ALLOC2D(Cinv  ,5,5); nr_identf(Cinv ,1,4,1,4);
  VIO_ALLOC2D(SR    ,5,5); nr_identf(SR   ,1,4,1,4);
  VIO_ALLOC2D(SRinv ,5,5); nr_identf(SRinv,1,4,1,4);
  VIO_ALLOC2D(Sinv  ,5,5); nr_identf(Sinv ,1,4,1,4); 
  VIO_ALLOC2D(T     ,5,5); nr_identf(T    ,1,4,1,4);
  VIO_ALLOC2D(Tinv  ,5,5); nr_identf(Tinv ,1,4,1,4);
  VIO_ALLOC2D(C     ,5,5); nr_identf(C    ,1,4,1,4);
  VIO_ALLOC2D(R     ,5,5); nr_identf(R    ,1,4,1,4);

  VIO_ALLOC2D(center_of_rotation ,5,5);        /* make column vectors */
  VIO_ALLOC2D(result             ,5,5);
  VIO_ALLOC2D(unit_vec           ,5,5);
  VIO_ALLOC2D(x                  ,5,5);
  VIO_ALLOC2D(y                  ,5,5);
  VIO_ALLOC2D(z                  ,5,5);
  VIO_ALLOC2D(nz                 ,5,5);
  VIO_ALLOC2D(y_on_z             ,5,5);
  VIO_ALLOC2D(ortho_y            ,5,5);

  ALLOC(tmp ,4);
  ALLOC(ang ,4);

  for(i=0; i<=3; i++)        /* copy the input matrix */
    for(j=0; j<=3; j++)
      xmat[i+1][j+1] = (float)Transform_elem(*trans,i,j);

  /* -------------DETERMINE THE TRANSLATION FIRST! ---------  */

                                /* see where the center of rotation is displaced... */

  FILL_NR_COLVEC( center_of_rotation, center[0], center[1], center[2] );


  invertmatrix(4, xmat, TMP1);        /* get inverse of the matrix */

  matrix_multiply( 4, 4, 1, xmat, center_of_rotation, result); /* was TMP! in place of xmat */

  SUB_NR_COLVEC( result, result, center_of_rotation );

  for(i=0; i<=2; i++) 
    translations[i] = result[i+1][1];

  /* -------------NOW GET THE SCALING VALUES! ----------------- */

  for(i=0; i<=2; i++) 
    tmp[i+1] = -translations[i];
  translation_to_homogeneous(3, tmp, Tinv); 

  for(i=0; i<=2; i++) 
    tmp[i+1] = center[i];
  translation_to_homogeneous(3, tmp, C); 
  for(i=0; i<=2; i++) 
    tmp[i+1] = -center[i];
  translation_to_homogeneous(3, tmp, Cinv); 


  matrix_multiply(4,4,4, xmat, C, TMP1);    /* get scaling*shear*rotation matrix */


  matrix_multiply(4,4,4, Tinv, TMP1, TMP1);


  matrix_multiply(4,4,4, Cinv, TMP1,    SR);

  invertmatrix(4, SR, SRinv);        /* get inverse of scaling*shear*rotation */

                                /* find each scale by mapping a unit vector backwards,
                                   and finding the magnitude of the result. */
  FILL_NR_COLVEC( unit_vec, 1.0, 0.0, 0.0 );
  matrix_multiply( 4, 4, 1, SRinv, unit_vec, result);
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
  matrix_multiply( 4, 4, 1, SRinv, unit_vec, result);
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
  matrix_multiply( 4, 4, 1, SRinv, unit_vec, result);
  magnitude = MAG_NR_COLVEC( result );
  if (magnitude != 0.0) {
    scales[2] = 1/magnitude;
    Sinv[3][3] = magnitude;
  }
  else {
    scales[2] = 1.0;
    Sinv[3][3] = 1.0;
  }

  /* ------------NOW GET THE SHEARS, using the info from above ----- */

  /* SR contains the [scale][shear][rot], must multiply [inv(scale)]*SR
     to get shear*rot. */

                                /* make [scale][rot] */
  matrix_multiply(4,4, 4, Sinv, SR,  TMP1);

                                /* get inverse of [scale][rot] */
  invertmatrix(4, TMP1, SRinv);        


  FILL_NR_COLVEC(x, SRinv[1][1], SRinv[2][1], SRinv[3][1]);
  FILL_NR_COLVEC(y, SRinv[1][2], SRinv[2][2], SRinv[3][2]);
  FILL_NR_COLVEC(z, SRinv[1][3], SRinv[2][3], SRinv[3][3]);



                                /* get normalized z direction  */
  magz = MAG_NR_COLVEC(z);
  SCALAR_MULT_NR_COLVEC( nz, z, 1/magz );

                                /* get a direction perpendicular 
                                   to z, in the yz plane.  */
  scalar = DOTSUM_NR_COLVEC(  y, nz );
  SCALAR_MULT_NR_COLVEC( y_on_z, nz, scalar );

  SUB_NR_COLVEC( result, y, y_on_z ); /* result = y - y_on_z */
  scalar = MAG_NR_COLVEC( result);     /* ortho_y = result ./ norm(result)  */
  SCALAR_MULT_NR_COLVEC( ortho_y, result, 1/scalar);


                                /* GET C for the skew matrix */

  scalar = DOTSUM_NR_COLVEC( y, nz ); /* project y onto z */
  magy   = MAG_NR_COLVEC(y);
  ci = scalar / sqrt((double)( magy*magy - scalar*scalar)) ;
                                /* GET B for the skew matrix */

                                /*    first need a1 */

  a1     = DOTSUM_NR_COLVEC( x, ortho_y ); /* project x onto ortho_y */
  magx   = MAG_NR_COLVEC(x);

                                /*    now get B  */

  scalar = DOTSUM_NR_COLVEC( x, nz );
  bi = scalar / sqrt((double)( magx*magx - scalar*scalar - a1*a1)) ;

                                /* GET A for skew matrix  */

  ai = a1 / sqrt((double)( magx*magx - scalar*scalar - a1*a1));

                                /* normalize the inverse shear parameters.
                                   so that there is no scaling in the matrix 
                                   (all scaling is already accounted for 
                                   in sx,sy and sz. */

  n1 = sqrt((double)(1 + ai*ai + bi*bi));
  n2 = sqrt((double)(1 + ci*ci));

  ai = ai / n1;
  bi = bi / n1;
  ci = ci / n2;

                                /* ai,bi,c1 now contain the parameters for 
                                   the inverse NORMALIZED shear matrix 
                                   (i.e., norm(col_i) = 1.0). */

  
  /* ------------NOW GET THE ROTATION ANGLES!----- */

                                  /*  since SR = [scale][shear][rot], then
                                    rot = [inv(shear)][inv(scale)][SR] */

  nr_identf(TMP1 ,1,4,1,4);        /* make inverse scale matrix */
  TMP1[1][1] = 1/scales[0];
  TMP1[2][2] = 1/scales[1];
  TMP1[3][3] = 1/scales[2];

  nr_identf(TMP2 ,1,4,1,4);        /* make_inverse normalized shear matrix */
  TMP2[1][1] = sqrt((double)(1 - ai*ai - bi*bi));
  TMP2[2][2] = sqrt((double)(1 - ci*ci));
  TMP2[2][1] = ai;
  TMP2[3][1] = bi;
  TMP2[3][2] = ci;


                                /* extract rotation matrix */
  matrix_multiply(4,4, 4, TMP2, TMP1, T);
  matrix_multiply(4,4, 4, T,    SR,   R);

                                /* get rotation angles */
  if (!rotmat_to_ang(R, ang)) {
    (void)fprintf(stderr,"Cannot convert R to radians!");
    printmatrix(3,3,R);
    return(FALSE);
  }

  for(i=0; i<=2; i++)
    rotations[i] = ang[i+1];

  /* ------------NOW ADJUST THE SCALE AND SKEW PARAMETERS ------------ */

  invertmatrix(4, T, Tinv);        /* get inverse of the matrix */
  

  scales[0] = Tinv[1][1];
  scales[1] = Tinv[2][2];
  scales[2] = Tinv[3][3];
  shears[0] = Tinv[2][1]/scales[1] ;
  shears[1] = Tinv[3][1]/scales[2] ;
  shears[2] = Tinv[3][2]/scales[2] ;
  

  VIO_FREE2D(xmat);
  VIO_FREE2D(TMP1);
  VIO_FREE2D(TMP2);
  VIO_FREE2D(Cinv);
  VIO_FREE2D(SR  );
  VIO_FREE2D(SRinv);
  VIO_FREE2D(Sinv);
  VIO_FREE2D(T   );
  VIO_FREE2D(Tinv);
  VIO_FREE2D(C   );
  VIO_FREE2D(R   );
  
  VIO_FREE2D(center_of_rotation);
  VIO_FREE2D(result            );
  VIO_FREE2D(unit_vec          );
  VIO_FREE2D(x                 );
  VIO_FREE2D(y                 );
  VIO_FREE2D(z                 );
  VIO_FREE2D(nz                );
  VIO_FREE2D(y_on_z            );
  VIO_FREE2D(ortho_y           );

  FREE(ang);
  FREE(tmp);

  return(TRUE);
}


/*
   function getparams will get the trnsform  parameters from a 
   transformation matrix 'tran' (that has already had the translation 
   componants removed).
   in description below, I assume that tran is a forward xform, from 
   native space to talairach space.
   I assume that trans is a 4x4 homogeneous matrix, with the 
   principal axis stored in the upper left 3x3 matrix.  The first col
   of tran represents is the coordinate of (1,0,0)' mapped through the
   transformation.  (this means  vec2 = tran * vec1).
  
   trans = [scale][shear][rot]
         = [scale][shear][rz][ry][rx];

   the shear matrix is constrained to be of the form:
     shear = [1 0 0 0
              f 1 0 0
              g h 1 0
              0 0 0 1];
     where f,g,h can take on any value.

   the scale matrix is constrained to be of the form:
     scale = [sx 0  0  0
              0  sy 0  0
              0  0  sz 0
              0  0  0  1];
     where sx,sy,sz must be positive.
   
  all rotations are assumed to be in the range -pi/2..pi/2

  the rotation angles are returned as radians and are applied 
  counter clockwise, when looking down the axis (from the positive end
  towards the origin).

  trans is assumed to be invertible.

  i assume a coordinate system:
             ^ y
             |
             |
             |
             |_______> x
            /
           /
          /z  (towards the viewer).

  
  
  procedure: 
          start with t = inv(tran) (after removing translation from tran )
  
          t = inv(r)*inv(sh)*inv(sc)
                  t maps the talairach space unit vectors into native
                space, which means that the columns of T are the
                direction cosines of these vectors x,y,z

   (1)  the length of the vectors x,y,z give the inverse scaling
        parameters:
             sx = 1/norm(x); sy = 1/norm(y); sz = 1/norm(z);

   (2)  with the constraints on the form of sh above, inv(sh) has the
        same form.  let inv(sh) be as above in terms of a,b and c.
          inv(sh) = [1 0 0 0; a 1 0 0; b c 1 0; 0 0 0 1];

        for c: project y onto z and normalize:

                  /     |y.z|^2     \(1/2)
             c = <  ---------------  >
                  \ |y|^2 - |y.z|^2 /

        for b: project x onto z and normalize

                  /        |x.z|^2        \(1/2)
             b = <  ---------------------  >
                  \ |x|^2 - |x.z|^2 - a^2 /
     
          where a is the projection of x onto the coordinate sys Y axis.

        for a: project x onto z and normalize
           a is taken from (b) above, and normalized... see below

  (3) rots are returned by getrots by giving the input transformation:
        rot_mat = [inv(sh)][inv(sc)][trans]

  (4) once completed, the parameters of sx,sy,sz and a,b,c are
      adjusted so that they maintain the matrix contraints above.


*/



VIO_BOOL extract2_parameters_from_matrix_quater(VIO_Transform *trans,
                                                      double *center,
                                                      double *translations,
                                                      double *scales,
                                                      double *shears,
                                                      double *quaternions)
{
  int 
    i,j;

  float 
    n1,n2,
    magnitude, magz, magx, magy, ai,bi,ci,scalar,a1,
    **center_of_rotation,
    **result,
    **unit_vec,
    *ang,*tmp,**x,**y,**z, **nz, **y_on_z, **ortho_y,
    **xmat,**T,**Tinv,**C,**Sinv,
    **R,**SR,**SRinv,**Cinv,**TMP1,**TMP2;

  VIO_ALLOC2D(xmat  ,5,5); nr_identf(xmat ,1,4,1,4);
  VIO_ALLOC2D(TMP1  ,5,5); nr_identf(TMP1 ,1,4,1,4);
  VIO_ALLOC2D(TMP2  ,5,5); nr_identf(TMP2 ,1,4,1,4);
  VIO_ALLOC2D(Cinv  ,5,5); nr_identf(Cinv ,1,4,1,4);
  VIO_ALLOC2D(SR    ,5,5); nr_identf(SR   ,1,4,1,4);
  VIO_ALLOC2D(SRinv ,5,5); nr_identf(SRinv,1,4,1,4);
  VIO_ALLOC2D(Sinv  ,5,5); nr_identf(Sinv ,1,4,1,4); 
  VIO_ALLOC2D(T     ,5,5); nr_identf(T    ,1,4,1,4);
  VIO_ALLOC2D(Tinv  ,5,5); nr_identf(Tinv ,1,4,1,4);
  VIO_ALLOC2D(C     ,5,5); nr_identf(C    ,1,4,1,4);
  VIO_ALLOC2D(R     ,5,5); nr_identf(R    ,1,4,1,4);

  VIO_ALLOC2D(center_of_rotation ,5,5);        /* make column vectors */
  VIO_ALLOC2D(result             ,5,5);
  VIO_ALLOC2D(unit_vec           ,5,5);
  VIO_ALLOC2D(x                  ,5,5);
  VIO_ALLOC2D(y                  ,5,5);
  VIO_ALLOC2D(z                  ,5,5);
  VIO_ALLOC2D(nz                 ,5,5);
  VIO_ALLOC2D(y_on_z             ,5,5);
  VIO_ALLOC2D(ortho_y            ,5,5);

  ALLOC(tmp ,4);
  ALLOC(ang ,4);

  for(i=0; i<=3; i++)        /* copy the input matrix */
    for(j=0; j<=3; j++)
      xmat[i+1][j+1] = (float)Transform_elem(*trans,i,j);
  
  /* -------------DETERMINE THE TRANSLATION FIRST! ---------  */

                                /* see where the center of rotation is displaced... */

  FILL_NR_COLVEC( center_of_rotation, center[0], center[1], center[2] );


  invertmatrix(4, xmat, TMP1);        /* get inverse of the matrix */

  matrix_multiply( 4, 4, 1, xmat, center_of_rotation, result); /* was TMP1 in place of xmat */

  SUB_NR_COLVEC( result, result, center_of_rotation );

  for(i=0; i<=2; i++) 
    translations[i] = result[i+1][1];

  /* -------------NOW GET THE SCALING VALUES! ----------------- */

  for(i=0; i<=2; i++) 
    tmp[i+1] = -translations[i];
  translation_to_homogeneous(3, tmp, Tinv); 

  for(i=0; i<=2; i++) 
    tmp[i+1] = center[i];
  translation_to_homogeneous(3, tmp, C); 
  for(i=0; i<=2; i++) 
    tmp[i+1] = -center[i];
  translation_to_homogeneous(3, tmp, Cinv); 


  matrix_multiply(4,4,4, xmat, C, TMP1);    /* get scaling*shear*rotation matrix */

 
  matrix_multiply(4,4,4, Tinv, TMP1, TMP1);


  matrix_multiply(4,4,4, Cinv, TMP1,    SR);

  invertmatrix(4, SR, SRinv);        /* get inverse of scaling*shear*rotation */
  
 
                                /* find each scale by mapping a unit vector backwards,
                                   and finding the magnitude of the result. */
  FILL_NR_COLVEC( unit_vec, 1.0, 0.0, 0.0 );
  matrix_multiply( 4, 4, 1, SRinv, unit_vec, result);
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
  matrix_multiply( 4, 4, 1, SRinv, unit_vec, result);
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
  matrix_multiply( 4, 4, 1, SRinv, unit_vec, result);
  magnitude = MAG_NR_COLVEC( result );
  
  if (magnitude != 0.0) {
    scales[2] = 1/magnitude;
    Sinv[3][3] = magnitude;
  }
  else {
    scales[2] = 1.0;
    Sinv[3][3] = 1.0;
  }

  /* ------------NOW GET THE SHEARS, using the info from above ----- */

  /* SR contains the [scale][shear][rot], must multiply [inv(scale)]*SR
     to get shear*rot. */

                                /* make [scale][rot] */
  matrix_multiply(4,4, 4, Sinv, SR,  TMP1);

                                /* get inverse of [scale][rot] */
  invertmatrix(4, TMP1, SRinv);        


  FILL_NR_COLVEC(x, SRinv[1][1], SRinv[2][1], SRinv[3][1]);
  FILL_NR_COLVEC(y, SRinv[1][2], SRinv[2][2], SRinv[3][2]);
  FILL_NR_COLVEC(z, SRinv[1][3], SRinv[2][3], SRinv[3][3]);



                                /* get normalized z direction  */
  magz = MAG_NR_COLVEC(z);
  SCALAR_MULT_NR_COLVEC( nz, z, 1/magz );

                                /* get a direction perpendicular 
                                   to z, in the yz plane.  */
  scalar = DOTSUM_NR_COLVEC(  y, nz );
  SCALAR_MULT_NR_COLVEC( y_on_z, nz, scalar );

  SUB_NR_COLVEC( result, y, y_on_z ); /* result = y - y_on_z */
  scalar = MAG_NR_COLVEC( result);     /* ortho_y = result ./ norm(result)  */
  SCALAR_MULT_NR_COLVEC( ortho_y, result, 1/scalar);


                                /* GET C for the skew matrix */

  scalar = DOTSUM_NR_COLVEC( y, nz ); /* project y onto z */
  magy   = MAG_NR_COLVEC(y);
  ci = scalar / sqrt((double)( magy*magy - scalar*scalar)) ;
                                /* GET B for the skew matrix */

                                /*    first need a1 */

  a1     = DOTSUM_NR_COLVEC( x, ortho_y ); /* project x onto ortho_y */
  magx   = MAG_NR_COLVEC(x);

                                /*    now get B  */

  scalar = DOTSUM_NR_COLVEC( x, nz );
  bi = scalar / sqrt((double)( magx*magx - scalar*scalar - a1*a1)) ;

                                /* GET A for skew matrix  */

  ai = a1 / sqrt((double)( magx*magx - scalar*scalar - a1*a1));

                                /* normalize the inverse shear parameters.
                                   so that there is no scaling in the matrix 
                                   (all scaling is already accounted for 
                                   in sx,sy and sz. */

  n1 = sqrt((double)(1 + ai*ai + bi*bi));
  n2 = sqrt((double)(1 + ci*ci));

  ai = ai / n1;
  bi = bi / n1;
  ci = ci / n2;

                                /* ai,bi,c1 now contain the parameters for 
                                   the inverse NORMALIZED shear matrix 
                                   (i.e., norm(col_i) = 1.0). */

  
  /* ------------NOW GET THE ROTATION ANGLES!----- */

                                  /*  since SR = [scale][shear][rot], then
                                    rot = [inv(shear)][inv(scale)][SR] */

  nr_identf(TMP1 ,1,4,1,4);        /* make inverse scale matrix */
  TMP1[1][1] = 1/scales[0];
  TMP1[2][2] = 1/scales[1];
  TMP1[3][3] = 1/scales[2];

  nr_identf(TMP2 ,1,4,1,4);        /* make_inverse normalized shear matrix */
  TMP2[1][1] = sqrt((double)(1 - ai*ai - bi*bi));
  TMP2[2][2] = sqrt((double)(1 - ci*ci));
  TMP2[2][1] = ai;
  TMP2[3][1] = bi;
  TMP2[3][2] = ci;


                                /* extract rotation matrix */
  matrix_multiply(4,4, 4, TMP2, TMP1, T);
  matrix_multiply(4,4, 4, T,    SR,   R);

 
                                /* get rotation angles */
  extract_quaternions(R,quaternions);

  /* ------------NOW ADJUST THE SCALE AND SKEW PARAMETERS ------------ */

  invertmatrix(4, T, Tinv);        /* get inverse of the matrix */
  

  scales[0] = Tinv[1][1];
  scales[1] = Tinv[2][2];
  scales[2] = Tinv[3][3];
  shears[0] = Tinv[2][1]/scales[1] ;
  shears[1] = Tinv[3][1]/scales[2] ;
  shears[2] = Tinv[3][2]/scales[2] ;
  

  VIO_FREE2D(xmat);
  VIO_FREE2D(TMP1);
  VIO_FREE2D(TMP2);
  VIO_FREE2D(Cinv);
  VIO_FREE2D(SR  );
  VIO_FREE2D(SRinv);
  VIO_FREE2D(Sinv);
  VIO_FREE2D(T   );
  VIO_FREE2D(Tinv);
  VIO_FREE2D(C   );
  VIO_FREE2D(R   );
  
  VIO_FREE2D(center_of_rotation);
  VIO_FREE2D(result            );
  VIO_FREE2D(unit_vec          );
  VIO_FREE2D(x                 );
  VIO_FREE2D(y                 );
  VIO_FREE2D(z                 );
  VIO_FREE2D(nz                );
  VIO_FREE2D(y_on_z            );
  VIO_FREE2D(ortho_y           );

  FREE(ang);
  FREE(tmp);

  return(TRUE);
}















