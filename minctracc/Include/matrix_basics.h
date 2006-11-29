/* ----------------------------- MNI Header -----------------------------------
@NAME       : matrix_basics.h
@DESCRIPTION: File containing declarations of routines for doing basic 
              matrix calculations
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

@CREATED    : January 31, 1992 (Peter Neelin)
@MODIFIED   : Mon Jun  7 11:41:02 EST 1993 (LC) update with new routines
---------------------------------------------------------------------------- */

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

