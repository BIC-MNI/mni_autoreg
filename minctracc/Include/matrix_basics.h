/* ----------------------------- MNI Header -----------------------------------
@NAME       : matrix_basics.h
@DESCRIPTION: File containing declarations of routines for doing basic 
              matrix calculations
@CREATED    : January 31, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */

public void printmatrix(int rows, int cols, float **the_matrix);
public void calc_centroid(int npoints, int ndim, float **points, 
                          float *centroid);
public void translate(int npoints, int ndim, float **points, 
                      float *translation, float **newpoints);
public void transpose(int rows, int cols, float **mat, float **mat_transpose);
public void matrix_multiply(int ldim, int mdim, int ndim, 
                            float **Amat, float **Bmat, float **Cmat);
public float trace(int size, float **the_matrix);
public void matrix_scalar_multiply(int rows, int cols, float scalar, 
                            float **the_matrix, float **product);

