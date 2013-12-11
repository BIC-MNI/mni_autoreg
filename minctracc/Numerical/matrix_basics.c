/* ----------------------------- MNI Header -----------------------------------
@NAME       : matrix_basics.c
@DESCRIPTION: File containing routines for doing basic matrix calculations
@METHOD     : Contains routines :
                 printmatrix
                 calc_centroid
                 translate
                 transpose
                 matrix_multiply
                 trace
                 matrix_scalar_multiply
                 invertmatrix
@CALLS      : 
@COPYRIGHT  :
              Copyright 1993 Peter Neelin, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.

@CREATED    : January 31, 1992 (Peter Neelin)
@MODIFIED   :  $Log: matrix_basics.c,v $
@MODIFIED   :  Revision 96.6  2006-11-29 09:09:33  rotor
@MODIFIED   :   * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   :  Revision 96.5  2005/07/20 20:45:49  rotor
@MODIFIED   :      * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :      * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :      * Removed old VOLUME_IO cruft #defines
@MODIFIED   :      * Fixed up all Makefile.am's in subdirs
@MODIFIED   :      * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :      * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   :  Revision 96.4  2004/02/12 05:54:27  rotor
@MODIFIED   :   * removed /static defs
@MODIFIED   :
@MODIFIED   :  Revision 96.3  2002/12/13 21:16:30  lenezet
@MODIFIED   :  nonlinear in 2D has changed. The option -2D-non-lin is no more necessary. The grid transform has been adapted to feet on the target volume whatever is size. The Optimization is done on the dimensions for which "count" is greater than 1.
@MODIFIED   :
@MODIFIED   :  Revision 96.2  2002/03/26 14:15:40  stever
@MODIFIED   :  Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   :  Revision 96.1  1997/11/03 19:59:49  louis
@MODIFIED   :  - now include volume_io/internal_volume_io.h instead of volume_io.h
@MODIFIED   :
 * Revision 96.0  1996/08/21  18:21:58  louis
 * Release of MNI_AutoReg version 0.96
 *
 * Revision 9.6  1996/08/21  18:21:54  louis
 * Pre-release
 *
 * Revision 9.5  1996/08/12  14:15:50  louis
 * Never released version 0.95
 *
 * Revision 1.10  1996/08/12  14:15:46  louis
 * Pre-release
 *
 * Revision 1.9  1995/09/11  12:37:16  collins
 * All refs to numerical recipes routines have been replaced.
 * this is an updated working version - corresponds to mni_reg-0.1g
 * \
 *
 * Revision 1.8  1995/02/22  08:56:06  collins
 * Montreal Neurological Institute version.
 * compiled and working on SGI.  this is before any changes for SPARC/
 * Solaris.
 *
 * Revision 1.7  94/04/06  11:48:41  louis
 * working linted version of linear + non-linear registration based on Lvv
 * operator working in 3D
 * 
 * Revision 1.6  94/02/21  16:35:44  louis
 * version before feb 22 changes
 * 
 * Revision 1.5  93/11/15  16:27:04  louis
 * working version, with new library, with RCS revision stuff,
 * before deformations included
 * 

Tue Jun  1 12:46:44 EST 1993 LC
     added routines to make identity matrices, copy matrices, make rotation matrices
     each routine has two copies (one for float, one for double parameters)

Fri Jun  4 14:10:34 EST 1993 LC
     moved the *to_homogeneous() routines here, from procrustes.c
     changed all homogeneous transformation calls, so that the 
       homogeneous transformation has the translations< on the right,
       ie new_vec = matrix * old_vec!

---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /static-cvsroot/registration/mni_autoreg/minctracc/Numerical/matrix_basics.c,v 96.6 2006-11-29 09:09:33 rotor Exp $";
#endif

#include <config.h>
#include <volume_io.h>
#include "local_macros.h"

/* external calls: */
void   make_rots(float **xmat, float data_rot_x, float data_rot_y, float data_rot_z);

/* Routines defined in this file */

void printmatrix(int rows, int cols, float **the_matrix);

void calc_centroid(int npoints, int ndim, float **points, 
                          float *centroid);

void translate(int npoints, int ndim, float **points, 
                      float *translation, float **newpoints);

void transpose(int rows, int cols, float **mat, float **mat_transpose);

void invertmatrix(int n, float **mat, float **mat_invert);

void raw_matrix_multiply(int ldim, int mdim, int ndim, 
                                float **Amat, float **Bmat, float **Cmat);

void matrix_multiply(int ldim, int mdim, int ndim, 
                            float **Amat, float **Bmat, float **Cmat);

float trace(int size, float **the_matrix);

void matrix_scalar_multiply(int rows, int cols, float scalar, 
                            float **the_matrix, float **product);

void nr_identd(double **A, int m1, int m2, int n1, int n2 );
void nr_identf(float **A, int m1, int m2, int n1, int n2 );

void nr_copyd(double **A, int m1, int m2, int n1, int n2, double **B );
void nr_copyf(float  **A, int m1, int m2, int n1, int n2, float **B );

void nr_rotxd(double **M, double a);
void nr_rotxf(float **M, float a);

void nr_rotyd(double **M,double a);
void nr_rotyf(float **M, float a);

void nr_rotzd(double **M,double a);
void nr_rotzf(float **M, float a);

void nr_multd(double **A, int mA1, int mA2, int nA1, int nA2,
                     double **B, int mB1, int mB2, int nB1, int nB2, 
                     double **C);
void nr_multf(float **A, int mA1, int mA2, int nA1, int nA2,
                     float **B, int mB1, int mB2, int nB1, int nB2, 
                     float **C);


void transformations_to_homogeneous(int ndim, 
                  float *translation, float *centre_of_rotation,
                  float **rotation, float scale,
                  float **transformation);

void translation_to_homogeneous(int ndim, float *translation,
                                       float **transformation);

void rotation_to_homogeneous(int ndim, float **rotation,
                                       float **transformation);



/* ----------------------------- MNI Header -----------------------------------
@NAME       : printmatrix
@INPUT      : rows   - number of rows in matrix
              cols   - number of columns in matrix
              the_matrix - matrix to be printed (in zero offset form).
                 The dimensions of this matrix should be defined to be 
                 1 to rows and 1 to cols.
@OUTPUT     : (nothing)
@RETURNS    : (nothing)
@DESCRIPTION: Prints out a matrix on stdout with one row per line.
@METHOD     : 
@GLOBALS    : (none)
@CALLS      : (nothing special)
@CREATED    : Feb. 26, 1990 (Weiqian Dai)
@MODIFIED   : January 31, 1992 (Peter Neelin)
                 - change to roughly NIL-abiding code
---------------------------------------------------------------------------- */
void printmatrix(int rows, int cols, float **the_matrix)
{
   int i,j;
   float f;

   /* Loop through rows and columns, printing one row per line */
   for (i=1; i <= rows; ++i) {
      for (j=1; j <= cols; ++j) {
         f=the_matrix[i][j];
         (void) print(" %10.6f ",f);
      }
      (void) print("\n");
   }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : calc_centroid
@INPUT      : npoints - number of points
              ndim    - number of dimensions
              points  - points matrix (in zero offset form).
                 The dimensions of this matrix should be defined to be 
                 1 to npoints and 1 to ndim.
@OUTPUT     : centroid - vector of centroid of points (in num. rec. form)
                 This vector should run from 1 to ndim.
@RETURNS    : (nothing)
@DESCRIPTION: Calculates the centroid of a number of points in ndim dimensions.
@METHOD     : 
@GLOBALS    : (none)
@CALLS      : (nothing special)
@CREATED    : Feb. 26, 1990 (Weiqian Dai)
@MODIFIED   : January 31, 1992 (Peter Neelin)
                 - change to roughly NIL-abiding code and modified calling
                 sequence.
---------------------------------------------------------------------------- */
void calc_centroid(int npoints, int ndim, float **points, 
                          float *centroid)
{
   int i,j;

   /* Loop over each dimension */
   for (i=1; i <= ndim; ++i) {
      /* Add up points and divide by number of points */
      centroid[i] = 0;
      for (j=1; j <= npoints; ++j) {
         centroid[i] += points[j][i];
      }
      if (npoints >0) centroid[i] /= (float) npoints;
   }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : translate
@INPUT      : npoints - number of points
              ndim    - number of dimensions
              points  - points matrix (in zero offset form).
                 The dimensions of this matrix should be defined to be 
                 1 to npoints and 1 to ndim.
              translation - translation vector (in num. rec. form, running
                 from 1 to ndim).
@OUTPUT     : newpoints - translated points matrix (see points). This matrix
                 can be the original points matrix.
@RETURNS    : (nothing)
@DESCRIPTION: Translates a set of points by a given translation.
@METHOD     : 
@GLOBALS    : (none)
@CALLS      : (nothing special)
@CREATED    : Feb. 26, 1990 (Weiqian Dai)
@MODIFIED   : January 31, 1992 (Peter Neelin)
                 - change to roughly NIL-abiding code and modified calling
                 sequence.
---------------------------------------------------------------------------- */
void translate(int npoints, int ndim, float **points, 
                      float *translation, float **newpoints)
{
   int i,j;

   for (i=1; i <= npoints; ++i) {
      for (j=1; j <= ndim; ++j) {
         newpoints[i][j] = points[i][j] + translation[j];
      }
   }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : transpose
@INPUT      : rows    - number of rows
              cols    - number of columns
              mat     - original matrix (in zero offset form).
                 The dimensions of this matrix should be defined to be 
                 1 to rows and 1 to cols.
@OUTPUT     : mat_transpose  - transposed matrix (in zero offset form,
                 with dimensions 1 to cols and 1 to rows). 
@RETURNS    : (nothing)
@DESCRIPTION: Transposes a matrix.
@METHOD     : 
@GLOBALS    : (none)
@CALLS      : (nothing special)
@CREATED    : Feb. 26, 1990 (Weiqian Dai)
@MODIFIED   : January 31, 1992 (Peter Neelin)
                 - change to roughly NIL-abiding code and modified calling
                 sequence.
Fri Jun  4 14:10:34 EST 1993 LC
    added the possibility to have input and out matrices the same!
---------------------------------------------------------------------------- */
void transpose(int rows, int cols, float **mat, float **mat_transpose)
{
   int i,j;

   float **Ctemp;

   if (mat==mat_transpose) {              /* if input and output the same, then alloc
                                         temporary space, so as not to overwrite
                                         the input before the compete transpose is 
                                         done. */
     /* Allocate a temporary matrix */
     VIO_ALLOC2D(Ctemp,cols+1,rows+1);
     
     for (i=1; i <= rows; ++i) {
       for (j=1; j <= cols; ++j) {
         Ctemp[j][i]=mat[i][j];
       }
     }
     
     /* Copy the result */
     for (i=1; i <= cols; ++i)
       for (j=1; j <= rows; ++j)
         mat_transpose[i][j] = Ctemp[i][j];
     
     /* Free the matrix */
     VIO_FREE2D(Ctemp);
   }
   else {
     for (i=1; i <= rows; ++i) {
       for (j=1; j <= cols; ++j) {
         mat_transpose[j][i]=mat[i][j];
       }
     }
   }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : invertmatrix
@INPUT      : n       - number of rows/or columns (must be square)
              mat     - original matrix (in zero offset form).
                 The dimensions of this matrix should be defined to be 
                 1 to n rows and 1 to n cols.
@OUTPUT     : mat_invert  - the inverted  matrix (in zero offset form,
                 with dimensions 1 to n cols and 1 to n rows). 
@RETURNS    : (nothing)
@DESCRIPTION: Inverts a matrix.
@METHOD     : 
@GLOBALS    : (none)
@CALLS      : (nothing special)
@CREATED    : Fri Jun  4 14:10:34 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */
void raw_invertmatrix(int n, float **mat, float **mat_invert)
{

  int 
    i,j;
  VIO_Real 
    **Rmat, **Rinv;

  VIO_ALLOC2D( Rmat, n, n );
  VIO_ALLOC2D( Rinv, n, n );

  for (i=1; i<=n; ++i)                /* copy the input matrix */
    for (j=1; j<=n; ++j) {
      Rmat[i-1][j-1] = mat[i][j];
    }

  (void)invert_square_matrix(n, Rmat, Rinv);

  for (i=1; i<=n; ++i)                /* copy the result */
    for (j=1; j<=n; ++j) {
      mat_invert[i][j] = Rinv[i-1][j-1];
    }

  VIO_FREE2D( Rmat );
  VIO_FREE2D( Rinv );

/*
                               this is the old inversion code
  float 
    d, **u, *col;
  int 
    i,j,*indx;


  u=mat rix(1,n,1,n);
  col=vec tor(1,n);
  indx=ivec tor(1,n);

  for (i=1; i<=n; ++i)                / * copy the input matrix * /
    for (j=1; j<=n; ++j)
      u[i][j] = mat[i][j];

  lud cmp(u,n,indx,&d);
  for(j=1; j<=n; ++j) {
    for(i=1; i<=n; ++i) col[i] = 0.0;
    col[j]=1.0;
    lub ksb(u,n,indx,col);
    for(i=1; i<=n; ++i) mat_invert[i][j]=col[i];
  }

  free_ matrix(u,1,n,1,n);
  free_ vector(col,1,n);
  free_ ivector(indx,1,n);

*/

}

void invertmatrix(int ndim, float **mat, float **mat_invert)
{
  float **Ctemp;
  int i,j;

  if (mat==mat_invert) {              /* if input and output the same, then alloc
                                         temporary space, so as not to overwrite
                                         the input as the inverse is being done. */
    /* Allocate a temporary matrix */
    VIO_ALLOC2D(Ctemp,ndim+1,ndim+1);
    
    /* invert the matrix */
    raw_invertmatrix(ndim, mat, Ctemp);
    
    /* Copy the result */
    for (i=1; i <= ndim; ++i)
      for (j=1; j <= ndim; ++j)
        mat_invert[i][j] = Ctemp[i][j];
    
    /* Free the matrix */
    VIO_FREE2D(Ctemp);
  }
  else {
    raw_invertmatrix(ndim, mat, mat_invert);
  }
}

/* ----------------------------- mni Header -----------------------------------
@NAME       : raw_matrix_multiply
@INPUT      : ldim, mdim, ndim - dimensions of matrices. Matrix Amat has
                 dimensions (ldim x mdim), matrix Bmat has dimension
                 (mdim x ndim) and resultant matrix has dimension
                 (ldim x ndim).
              Amat - First matrix of multiply (in zero offset form).
                 Dimensions are 1 to ldim and 1 to mdim.
              Bmat - Second matrix of multiply (in zero offset form).
                 Dimensions are 1 to mdim and 1 to ndim.
@OUTPUT     : Cmat - Resulting matrix (in zero offset form).
                 Dimensions are 1 to ldim and 1 to ndim. This matrix cannot
                 be either Amat or Bmat.
@RETURNS    : (nothing)
@DESCRIPTION: Multiplies two matrices.
@METHOD     : 
@GLOBALS    : (none)
@CALLS      : (nothing special)
@CREATED    : Feb. 26, 1990 (Weiqian Dai)
@MODIFIED   : January 31, 1992 (Peter Neelin)
                 - change to roughly NIL-abiding code and modified calling
                 sequence.
---------------------------------------------------------------------------- */
void raw_matrix_multiply(int ldim, int mdim, int ndim, 
                                float **Amat, float **Bmat, float **Cmat)
{
   int i,j,k;

   /* Zero the output matrix */
   for (i=1; i <= ldim; ++i)
      for (j=1; j <= ndim; ++j)
         Cmat[i][j]=0.;

   /* Calculate the product */
   for (i=1; i <= ldim; ++i)
      for (j=1; j <= ndim; ++j)
         for (k=1; k <=mdim; ++k)
            Cmat[i][j] += (Amat[i][k] * Bmat[k][j]);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : matrix_multiply
@INPUT      : ldim, mdim, ndim - dimensions of matrices. Matrix Amat has
                 dimensions (ldim x mdim), matrix Bmat has dimension
                 (mdim x ndim) and resultant matrix has dimension
                 (ldim x ndim).
              Amat - First matrix of multiply (in zero offset form).
                 Dimensions are 1 to ldim and 1 to mdim.
              Bmat - Second matrix of multiply (in zero offset form).
                 Dimensions are 1 to mdim and 1 to ndim.
@OUTPUT     : Cmat - Resulting matrix (in zero offset form).
                 Dimensions are 1 to ldim and 1 to ndim. This matrix can
                 be either matrix Amat or Bmat.
@RETURNS    : (nothing)
@DESCRIPTION: Multiplies two matrices.
@METHOD     : 
@GLOBALS    : (none)
@CALLS      : (nothing special)
@CREATED    : March 2, 1992 (Peter Neelin)
@MODIFIED   : March 2, 1992 (P.N.)
                 - Changed so that calling program can use an input matrix for
                 output.
---------------------------------------------------------------------------- */
void matrix_multiply(int ldim, int mdim, int ndim, 
                            float **Amat, float **Bmat, float **Cmat)
{
   int i,j;
   float **Ctemp;

   /* Allocate a temporary matrix */
   VIO_ALLOC2D(Ctemp, ldim+1, ndim+1);

   /* Do the multiplication */
   raw_matrix_multiply(ldim,mdim,ndim,Amat,Bmat,Ctemp);

   /* Copy the result */
   for (i=1; i <= ldim; ++i)
      for (j=1; j <= ndim; ++j)
         Cmat[i][j] = Ctemp[i][j];

   /* Free the matrix */
   VIO_FREE2D(Ctemp);
}
                  

/* ----------------------------- MNI Header -----------------------------------
@NAME       : trace
@INPUT      : size   - size of the_matrix (the_matrix should be square)
              the_matrix - matrix for which trace should be calculated (in 
                 zero offset form). Dimensions are 1 to size and 
                 1 to size.
@OUTPUT     : (none)
@RETURNS    : trace of matrix
@DESCRIPTION: Calculates the trace of a matrix.
@METHOD     : 
@GLOBALS    : (none)
@CALLS      : (nothing special)
@CREATED    : Feb. 26, 1990 (Weiqian Dai)
@MODIFIED   : January 31, 1992 (Peter Neelin)
                 - change to roughly NIL-abiding code and modified calling
                 sequence.
---------------------------------------------------------------------------- */
float trace(int size, float **the_matrix)
{
   float sum=0.;
   int i;

   for (i=1; i <= size; ++i) {
      sum += the_matrix[i][i];
   }

   return(sum);
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : matrix_scalar_multiply
@INPUT      : rows    - number of rows of the_matrix.
              cols    - number of columns of the_matrix
              scalar  - scalar by which the_matrix should be multiplied.
              the_matrix  - matrix to be multiplied (in zero offset 
                 form). Dimensions are 1 to rows and 1 to cols.
@OUTPUT     : product - result of multiply ( in zero offset form).
                 Dimensions are 1 to rows and 1 to cols. This matrix
                 can be the input matrix.
@RETURNS    : (nothing)
@DESCRIPTION: Multiplies a matrix by a scalar.
@METHOD     : 
@GLOBALS    : (none)
@CALLS      : (nothing special)
@CREATED    : Feb. 26, 1990 (Weiqian Dai)
@MODIFIED   : January 31, 1992 (Peter Neelin)
                 - change to roughly NIL-abiding code and modified calling
                 sequence.
---------------------------------------------------------------------------- */
void matrix_scalar_multiply(int rows, int cols, float scalar, 
                            float **the_matrix, float **product)
{
   int i,j;

   for (i=1; i <= rows; ++i)
      for (j=1; j<=cols; ++j)
         product[i][j]=scalar*the_matrix[i][j];
}




/* ----------------------------- MNI Header -----------------------------------
@NAME       : nr_identd, nr_identf - make identity matrix
@INPUT      : A - pointer to matrix
              m1,m2 - row limits
              n1,n2 - col limits
              (matrix in zero offset form, allocated by calling routine)
@OUTPUT     : identiy matrix in A
@RETURNS    : (nothing)
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : (none)
@CALLS      : (nothing special)
@CREATED    : Tue Jun  1 12:49:21 EST 1993 (Louis Collins)
@MODIFIED   : 

---------------------------------------------------------------------------- */
void nr_identd(double **A, int m1, int m2, int n1, int n2 )
{

   int i,j;

   for (i=m1; i<=m2; ++i)
      for (j=n1; j<=n2; ++j) {
         if (i==j) 
            A[i][j] = 1.0;
         else
            A[i][j] = 0.0;
      }
   
}

void nr_identf(float **A, int m1, int m2, int n1, int n2 )
{

   int i,j;

   for (i=m1; i<=m2; ++i)
      for (j=n1; j<=n2; ++j) {
         if (i==j) 
            A[i][j] = 1.0;
         else
            A[i][j] = 0.0;
      }
   
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : nr_copyd, nr_copyf - copy matrix
@INPUT      : A - source matrix
              m1,m2 - row limits
              n1,n2 - col limits
              (matrix in zero offset form, allocated by calling routine)
@OUTPUT     : B - copy of A
@RETURNS    : (nothing)
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : (none)
@CALLS      : (nothing special)
@CREATED    : Tue Jun  1 12:49:21 EST 1993 (Louis Collins)
@MODIFIED   : 

---------------------------------------------------------------------------- */
void nr_copyd(double **A, int m1, int m2, int n1, int n2, double **B )
{
   int i,j;

   for (i=m1; i<=m2; ++i)
      for (j=n1; j<=n2; ++j)
         B[i][j] = A[i][j];
}

void nr_copyf(float  **A, int m1, int m2, int n1, int n2, float **B )
{
   int i,j;

   for (i=m1; i<=m2; ++i)
      for (j=n1; j<=n2; ++j)
         B[i][j] = A[i][j];
}



/* ----------------------------- MNI Header -----------------------------------
@NAME       : nr_rotxd,nr_rotxf - make rot X matrix
@INPUT      : M - 4x4 matrix
              a - rotation angle in radians
              (matrix in zero offset form, allocated by calling routine)
@OUTPUT     : modified matrix M
@RETURNS    : (nothing)
@DESCRIPTION: 
@METHOD     : 
   rx =[1   0      0      0 
        0  cos(a)  -sin(a) 0
        0 sin(a)  cos(a) 0
        0   0      0      1];
@GLOBALS    : (none)
@CALLS      : (nothing special)
@CREATED    : Tue Jun  1 12:49:21 EST 1993 (Louis Collins)
@MODIFIED   : Tue Jun  8 08:44:59 EST 1993 (LC) changed to mat*vec format
---------------------------------------------------------------------------- */
void nr_rotxd(double **M, double a)
{
   nr_identd(M,1,4,1,4);

   M[2][2] = cos(a);   M[2][3] = -sin(a);
   M[3][2] = sin(a);   M[3][3] = cos(a);
}


void nr_rotxf(float **M, float a)
{
   nr_identf(M,1,4,1,4);

   M[2][2] = cos((double)a);    M[2][3] = -sin((double)a);
   M[3][2] = sin((double)a);   M[3][3] = cos((double)a);
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : nr_rotyd,nr_rotyf - make rot Y matrix
@INPUT      : M - 4x4 matrix
              a - rotation angle in radians
              (matrix in zero offset form, allocated by calling routine)
@RETURNS    : (nothing)
@DESCRIPTION: 
@METHOD     : 
ry = [  cos(a)   0 sin(a)  0 
        0       1   0       0
        -sin(a)  0  cos(a)   0
        0   0      0      1];
@GLOBALS    : (none)
@CALLS      : (nothing special)
@CREATED    : Tue Jun  1 12:49:21 EST 1993 (Louis Collins)
@MODIFIED   : Tue Jun  8 08:44:59 EST 1993 (LC) changed to mat*vec format
---------------------------------------------------------------------------- */
void nr_rotyd(double **M,double a)
{

   nr_identd(M,1,4,1,4);

   M[1][1] = cos(a);   M[1][3] = sin(a);
   M[3][1] = -sin(a);   M[3][3] = cos(a);
}

void nr_rotyf(float **M, float a)
{

   nr_identf(M,1,4,1,4);

   M[1][1] = cos((double)a);   M[1][3] = sin((double)a);
   M[3][1] = -sin((double)a);   M[3][3] = cos((double)a);
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : nr_rotzd, nr_rotzf - make rot Z matrix
@INPUT      : M - 4x4 matrix
              a - rotation angle in radians
              (matrix in zero offset form, allocated by calling routine)
@RETURNS    : (nothing)
@DESCRIPTION: 
@METHOD     : 
rz = [cos(a)  -sin(a) 0  0
      sin(a) cos(a) 0  0
        0     0      1  0
        0     0      0  1];
@GLOBALS    : (none)
@CALLS      : (nothing special)
@CREATED    : Tue Jun  1 12:49:21 EST 1993 (Louis Collins)
@MODIFIED   : Tue Jun  8 08:44:59 EST 1993 (LC) changed to mat*vec format
---------------------------------------------------------------------------- */
void nr_rotzd(double **M,double a)
{

   nr_identd(M,1,4,1,4);

   M[1][1] = cos(a);   M[1][2] = -sin(a);
   M[2][1] = sin(a);  M[2][2] = cos(a);
}

void nr_rotzf(float **M, float a)
{

   nr_identf(M,1,4,1,4);

   M[1][1] = cos((double)a);   M[1][2] = -sin((double)a);
   M[2][1] = sin((double)a);  M[2][2] = cos((double)a);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : nr_multd, nr_multf - mult matrix
@INPUT      : A - source matrix
              mA1,mA2 - row limits of A
              nA1,nA2 - col limits of A
              B - source matrix
              mB1,mB2 - row limits of B
              nB1,nB2 - col limits of B
              (matrix in zero offset form, allocated by calling routine)
@OUTPUT     : C = A * B
@RETURNS    : (nothing)
@DESCRIPTION: 
   Routine multiplies matrices A*B to give C. A is a mA x nA matrix and
   B is a mB x nB matrix. The result is returned in C which is mA x nB.
@METHOD     : 
@GLOBALS    : (none)
@CALLS      : (nothing special)
@CREATED    : Tue Jun  1 12:49:21 EST 1993 (Louis Collins)
@MODIFIED   : 

---------------------------------------------------------------------------- */

void nr_multd(double **A, int mA1, int mA2, int nA1, int nA2, 
         double **B, int mB1, int mB2, int nB1, int nB2, 
         double **C )
{
   int i, j, k;

   for ( k = mA1; k <= mA2; k++ ) {
      for ( i = nB1; i <= nB2; i++ ) {
         C[k][i] = 0.0;
         for ( j = mB1; j <= mB2; j++ ) {
            C[k][i] += A[k][j] * B[j][i];
         }
      }
   }

   return;
}


void nr_multf(float **A, int mA1, int mA2, int nA1, int nA2, 
         float **B, int mB1, int mB2, int nB1, int nB2, 
         float **C)
{
   int i, j, k;

   for ( k = mA1; k <= mA2; k++ ) {
      for ( i = nB1; i <= nB2; i++ ) {
         C[k][i] = 0.0;
         for ( j = mB1; j <= mB2; j++ ) {
            C[k][i] += A[k][j] * B[j][i];
         }
      }
   }

   return;
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : transformations_to_homogeneous
@INPUT      : ndim    - number of dimensions
              translation - zero offset vector (1 to ndim) that 
                 specifies the translation to be applied first.
              centre_of_rotation - zero offset vector (1 to ndim) that
                 specifies the centre of rotation and scaling.
              rotation - zero offset matrix (1 to ndim by 1 to ndim) 
                 for rotation about centre_of_rotation (applied after 
                 translation). Note that this matrix need not only specify
                 rotation/reflexion - any ndim x ndim matrix will work.
              scale - Scalar value giving global scaling to be applied after
                 translation and rotation.
@OUTPUT     : transformation - zero offset matrix (1 to ndim+1 by
                 1 to ndim+1) specifying the transformation for homogeneous 
                 coordinates. To apply this transformation, a point
                 vector should be pre-multiplied by this matrix, with the
                 last coordinate of the ndim+1 point vector having value
                 one. The calling routine must allocate space for this
                 matrix.
@RETURNS    : (nothing)
@DESCRIPTION: Computes a transformation matrix in homogeneous coordinates
              given a translation, a rotation matrix (or other 
              non-homogeneous matrix) and a global scaling factor.
              Transformations are applied in that order.
@METHOD     : Apply the following operations (multiply from left to right):
                 1) Translate by translation
                 2) Translate by -centre_of_rotation
                 3) Rotate
                 4) Scale
                 5) Translate by centre_of_rotation
@GLOBALS    : (none)
@CALLS      : translation_to_homogeneous
              matrix_multiply
              matrix_scalar_multiply
@CREATED    : February 7, 1992 (Peter Neelin)
@MODIFIED   : 
Fri Jun  4 14:10:34 EST 1993  LC
   changed matrices, so that they must be applied by pre-multiplication:
      ie newvec = matrix * oldvec
---------------------------------------------------------------------------- */
void transformations_to_homogeneous(int ndim, 
                  float *translation, float *centre_of_rotation,
                  float **rotation, float scale,
                  float **transformation)
{
   int i;
   int size;
   float *centre_translate;
   float **trans1, **trans2;
   float **trans_temp, **rotation_and_scale;

   size=ndim+1;

   /* Allocate matrices and vectors */
   ALLOC(centre_translate,ndim+1);
   VIO_ALLOC2D(trans1 ,size+1, size+1);
   VIO_ALLOC2D(trans2 ,size+1, size+1);
   VIO_ALLOC2D(trans_temp,size+1, size+1); 
   VIO_ALLOC2D(rotation_and_scale,ndim+1, ndim+1);


   /* Construct translation matrix */
   translation_to_homogeneous(ndim, translation, trans1);


   /* Construct translation matrix for centre of rotation and
      apply it */
   for (i=1; i<=ndim; i++) centre_translate[i] = -centre_of_rotation[i];
   translation_to_homogeneous(ndim, centre_translate, trans_temp);
   matrix_multiply(size, size, size, trans1, trans_temp, trans2);


   /* Scale rotation matrix, then convert it to homogeneous coordinates and
      apply it */
   matrix_scalar_multiply(ndim, ndim, scale, rotation, rotation_and_scale);
   rotation_to_homogeneous(ndim, rotation_and_scale, trans_temp);
   matrix_multiply(size, size, size, trans2, trans_temp, trans1);


   /* Return to centre of rotation */
   translation_to_homogeneous(ndim, centre_of_rotation, trans_temp);
   matrix_multiply(size, size, size, trans1, trans_temp, transformation);


   /* Free matrices */
   FREE(  centre_translate);
   VIO_FREE2D(trans1);
   VIO_FREE2D(trans2);
   VIO_FREE2D(trans_temp);
   VIO_FREE2D(rotation_and_scale);

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : translation_to_homogeneous
@INPUT      : ndim    - number of dimensions
              translation - zero offset vector (1 to ndim) that 
                 specifies the translation.
@OUTPUT     : transformation - zero offset matrix (1 to ndim+1 by
                 1 to ndim+1) specifying the transformation for homogeneous 
                 coordinates. To apply this transformation, a point
                 vector should be pre-multiplied by this matrix, with the
                 last coordinate of the ndim+1 point vector having value
                 one. The calling routine must allocate space for this
                 matrix.
@RETURNS    : (nothing)
@DESCRIPTION: Computes a transformation matrix in homogeneous coordinates
              given a translation.
@METHOD     : 
@GLOBALS    : (none)
@CALLS      : 
@CREATED    : February 7, 1992 (Peter Neelin)
@MODIFIED   : 
Fri Jun  4 14:10:34 EST 1993  LC
   changed matrices, so that they must be applied by pre-multiplication:
      ie newvec = matrix * oldvec
---------------------------------------------------------------------------- */
void translation_to_homogeneous(int ndim, float *translation,
                                       float **transformation)
{
   int i,j;
   int size;

   size=ndim+1;

   /* Construct translation matrix */
   for (i=1; i<=ndim; i++) {
      for (j=1; j<=size; j++) {
         if (i == j) {
            transformation[i][j] = 1.0;
         }
         else {
            transformation[i][j] = 0.0;
         }
      }
   }
   for (j=1; j<=ndim; j++) {
      transformation[j][size] = translation[j];
   }

   transformation[size][size] = 1.0;

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : rotation_to_homogeneous
@INPUT      : ndim    - number of dimensions
              rotation - zero offset matrix (1 to ndim by 1 to ndim) 
                 for rotation about origin. Note that this matrix need not 
                 only specify rotation/reflexion - any ndim x ndim matrix 
                 will work.
@OUTPUT     : transformation - zero offset matrix (1 to ndim+1 by
                 1 to ndim+1) specifying the transformation for homogeneous 
                 coordinates. To apply this transformation, a point
                 vector should be pre-multiplied by this matrix, with the
                 last coordinate of the ndim+1 point vector having value
                 one. The calling routine must allocate space for this
                 matrix.
@RETURNS    : (nothing)
@DESCRIPTION: Computes a transformation matrix in homogeneous coordinates
              given a rotation matrix.
@METHOD     : 
@GLOBALS    : (none)
@CALLS      : 
@CREATED    : February 7, 1992 (Peter Neelin)
@MODIFIED   : 
Fri Jun  4 14:10:34 EST 1993  LC
   changed matrices, so that they must be applied by pre-multiplication:
      ie newvec = matrix * oldvec
---------------------------------------------------------------------------- */
void rotation_to_homogeneous(int ndim, float **rotation,
                                       float **transformation)
{
   int i,j;
   int size;

   size=ndim+1;

   /* Construct  matrix */
   for (i=1; i<=size; i++) {
      for (j=1; j<=size; j++) {
         if ((i==size) || (j==size)) {
            transformation[i][j] = 0.0;
         }
         else {
            transformation[i][j] = rotation[i][j];
         }
      }
   }

   transformation[size][size] = 1.0;

}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : angles_to_homogeneous
@INPUT      : ndim    - number of dimensions
              angles - zero offset array (1 to ndim)
                 for rotation angles (in radians) about origin. 
@OUTPUT     : transformation - zero offset matrix (1 to ndim+1 by
                 1 to ndim+1) specifying the transformation for homogeneous 
                 coordinates. To apply this transformation, a point
                 vector should be pre-multiplied by this matrix, with the
                 last coordinate of the ndim+1 point vector having value
                 one. The calling routine must allocate space for this
                 matrix.
@RETURNS    : (nothing)
@DESCRIPTION: Computes a transformation matrix in homogeneous coordinates
              given a rotation matrix.
@METHOD     : 
@GLOBALS    : (none)
@CALLS      : 
@CREATED    : Fri Jun  4 14:10:34 EST 1993  LC
@MODIFIED   : 
---------------------------------------------------------------------------- */
void angles_to_homogeneous(int ndim, float *angles,
                                  float **transformation)
{
   int i,j;
   int size;
   float **rot_matrix;


   size=ndim+1;

   VIO_ALLOC2D(rot_matrix,5,5);


   if (ndim==2 || ndim==3) {

     if (ndim==2)
       nr_rotzf(rot_matrix,*angles );
     else
       make_rots(rot_matrix, 
                 (float)(angles[0]),
                 (float)(angles[1]),
                 (float)(angles[2]));

     /* Construct  matrix */
     for (i=1; i<=size; i++) {
       for (j=1; j<=size; j++) {
         if ((i==size) || (j==size)) {
           transformation[i][j] = 0.0;
         }
         else {
           transformation[i][j] = rot_matrix[i][j];
         }
       }
     }
     transformation[size][size] = 1.0;

   }
   else {
     (void)fprintf (stderr,"Can't handle %d dimensions in angles_to_homogeneous()\n",ndim);
     (void)fprintf (stderr,"Error in %s, line %d\n",__FILE__,__LINE__);
     exit(-1);
   }


   VIO_FREE2D(rot_matrix);
}
